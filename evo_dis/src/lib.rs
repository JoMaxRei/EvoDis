use bevy::prelude::*;
use foodweb::Foodweb;
use position::Position;
use rand::{Rng, SeedableRng, seq::IteratorRandom};
use rand_chacha::ChaCha8Rng;
use settings::{SimulationSettings, Statistics};
use species::{Individual, Source, Species};

mod dispersal;
mod foodweb;
mod position;
mod settings;
mod speciation;
mod species;
mod survival;

pub fn plugin(app: &mut App) {
    let sampler = ChaCha8Rng::seed_from_u64(0);
    app.register_type::<Speciation>()
        .insert_resource(Sampler(sampler))
        .init_resource::<SimulationSettings>()
        .init_resource::<Statistics>()
        .insert_state(SimulationState::Main)
        .insert_resource(SimulationTime(0.0))
        .register_type::<SimulationSettings>()
        .register_type::<Species>()
        .register_type::<Individual>()
        .add_plugins((speciation::plugin, dispersal::plugin, survival::plugin))
        .add_systems(Startup, setup)
        .add_systems(Update, log)
        .add_systems(
            Update,
            (handle_event.run_if(in_state(SimulationState::Main)),),
        )
        .add_systems(Update, check_input.run_if(in_state(SimulationState::Wait)))
        .add_systems(OnEnter(SimulationState::Done), print_statistics);
}

#[derive(Resource)]
struct SimulationTime(f64);

#[derive(Clone, Debug, Eq, Hash, PartialEq, States)]
enum SimulationState {
    Main,
    Prey,
    Speciation,
    CheckSurvivors,
    Dispersal,
    Wait,
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
    let _ = info_span!("setup", name = "setup").entered();
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
                .spawn((
                    Name::from(format!("Foodweb {}-{}", x, y)),
                    Foodweb,
                    Position::new(x, y),
                ))
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
    let _ = info_span!("handle_event", name = "handle_event").entered();
    if time.0 >= settings.speciations_per_patch {
        info!("time is over");
        next_state.set(SimulationState::Done);
        return;
    }

    let mut total_dispersal_rate = 0;
    for individuum in &individuums {
        total_dispersal_rate += species.get(individuum.0).unwrap().dispersal_rate;
    }
    trace!("total dispersal rate {}", total_dispersal_rate);

    let total_speciation_rate =
        settings.speciation_rate_per_individuum * individuums.iter().count() as u64;
    trace!("total speciation rate {}", total_speciation_rate);

    time.0 += 1.0
        / (1.0 + total_dispersal_rate as f64 / total_speciation_rate as f64)
        / settings.grid_size.x as f64
        / settings.grid_size.y as f64;
    trace!("time is now {}", time.0);

    let event_is_speciation = (total_dispersal_rate + total_speciation_rate)
        * sampler.0.gen_range(0..=1)
        >= total_dispersal_rate;

    if event_is_speciation {
        next_state.set(SimulationState::Speciation);
    } else {
        next_state.set(SimulationState::Dispersal);
    }
}

fn log(mut events: EventReader<StateTransitionEvent<SimulationState>>) {
    for event in events.read() {
        debug!("Moving from {:?} to {:?}", event.exited, event.entered);
    }
}

fn check_input(
    // input: Res<ButtonInput<KeyCode>>,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    // if input.just_pressed(KeyCode::Space) {
    next_state.set(SimulationState::Main);
    // }
}

fn print_statistics(statistics: Res<Statistics>) {
    info!("Failed dispersals {}", statistics.failed_dispersals);
    info!("Dispersals {}", statistics.dispersals);
    info!("Speciations {}", statistics.speciations);
    info!(
        "Events {}",
        statistics.dispersals
            + statistics.speciations
            + statistics.died_individuals
            + statistics.died_species
    );
    info!(
        "failrate {}",
        statistics.failed_dispersals as f32 / statistics.dispersals as f32
    );
}
