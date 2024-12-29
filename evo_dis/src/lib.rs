use bevy::prelude::*;
use foodweb::Foodweb;
use rand::{Rng, SeedableRng, seq::IteratorRandom};
use rand_chacha::ChaCha8Rng;
use settings::SimulationSettings;
use species::{Individual, Source, Species};

mod foodweb;
mod settings;
mod speciation;
mod species;

pub fn plugin(app: &mut App) {
    let sampler = ChaCha8Rng::seed_from_u64(0);
    app.register_type::<Speciation>()
        .insert_resource(Sampler(sampler))
        .init_resource::<SimulationSettings>()
        .insert_state(SimulationState::Main)
        .insert_resource(SimulationTime(0.0))
        .register_type::<SimulationSettings>()
        .register_type::<Species>()
        .register_type::<Individual>()
        .add_plugins(speciation::plugin)
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                handle_event.run_if(in_state(SimulationState::Main)),
                handle_dispersal.run_if(in_state(SimulationState::Dispersal)),
            ),
        );
}

#[derive(Resource)]
struct SimulationTime(f64);

#[derive(Clone, Debug, Eq, Hash, PartialEq, States)]
enum SimulationState {
    Main,
    Speciation,
    CheckSurvivors,
    Dispersal,
    Done,
}

#[derive(Resource)]
struct Sampler(ChaCha8Rng);

#[derive(Resource, Reflect)]
#[reflect(Resource)]
struct Speciation {
    foodweb: Entity,
    species: Entity,
}

fn setup(settings: Res<SimulationSettings>, mut sampler: ResMut<Sampler>, mut commands: Commands) {
    let resource = commands
        .spawn((
            Species {
                bodymass: 0.0,
                feeding_center: 0.0,
                feeding_range: 0.0,
                first_occurence: 0.0,
                predator_strength: 1.0,
                dispersal_rate: 0,
            },
            Source,
            Name::from("Species Resource"),
        ))
        .id();

    let first_species = commands
        .spawn((
            Species {
                bodymass: settings.mean_bodymass_ratio_predator_prey,
                feeding_center: 0.0,
                feeding_range: sampler
                    .0
                    .gen_range(settings.feeding_range.min..=settings.feeding_range.max),
                first_occurence: 0.0,
                predator_strength: Species::calculate_predator_strength(
                    settings.speciation_rate_per_individuum,
                    settings.zero_crossing,
                    settings.initial_dispersal_rate,
                ),
                dispersal_rate: 1 * settings.speciation_rate_per_individuum,
            },
            Name::from("Species 0"),
        ))
        .id();

    for x in 0..settings.grid_size.x {
        for y in 0..settings.grid_size.y {
            commands
                .spawn((Name::from(format!("Foodweb {}-{}", x, y)), Foodweb))
                .with_children(|foodweb| {
                    foodweb.spawn((Individual(resource), Source));
                    foodweb.spawn(Individual(first_species));
                });
        }
    }
}

fn handle_event(
    mut time: ResMut<SimulationTime>,
    settings: Res<SimulationSettings>,
    mut sampler: ResMut<Sampler>,
    individuums: Query<&Individual, Without<Source>>,
    species: Query<&Species>,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    if time.0 >= 10.0 {
        next_state.set(SimulationState::Done);
        return;
    }

    let mut total_dispersal_rate = 0;
    for individuum in &individuums {
        total_dispersal_rate += species.get(individuum.0).unwrap().dispersal_rate;
    }
    debug!("total dispersal rate {}", total_dispersal_rate);

    let total_speciation_rate =
        settings.speciation_rate_per_individuum * individuums.iter().count() as u64;
    debug!("total speciation rate {}", total_speciation_rate);

    time.0 += 1.0
        / (1.0 + total_dispersal_rate as f64 / total_speciation_rate as f64)
        / settings.grid_size.x as f64
        / settings.grid_size.y as f64;
    info!("time is now {}", time.0);

    let event_is_speciation = (total_dispersal_rate + total_speciation_rate)
        * sampler.0.gen_range(0..=1)
        >= total_dispersal_rate;

    if event_is_speciation {
        next_state.set(SimulationState::Speciation);
    } else {
        next_state.set(SimulationState::Dispersal);
    }
}

fn handle_dispersal(mut next_state: ResMut<NextState<SimulationState>>) {
    next_state.set(SimulationState::Main);
}
