use bevy::prelude::*;
use foodweb::Foodweb;
use rand::{Rng, SeedableRng, seq::IteratorRandom};
use rand_chacha::ChaCha8Rng;
use settings::SimulationSettings;
use species::{Individual, Source, Species};

mod foodweb;
mod settings;
mod species;

pub fn plugin(app: &mut App) {
    let sampler = ChaCha8Rng::seed_from_u64(0);
    app.add_event::<EventType>().insert_resource(Sampler(sampler))
        .init_resource::<SimulationSettings>()
        .insert_state(SimulationState::Main)
        .insert_resource(SimulationTime(0.0))
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                handle_event.run_if(in_state(SimulationState::Main)),
                (handle_speciation, handle_new_species).chain().run_if(in_state(SimulationState::Speciation)),
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
    Dispersal,
}

#[derive(Resource)]
struct Sampler(ChaCha8Rng);

#[derive(Event)]
enum EventType {
    /// Gets fired when a new species is inserted into the command queue.
    /// .0 contains the foodweb where it was from, .1 contains the species itself
    Speciation(Entity, Entity),
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
                dispersal_rate: 1,
            },
            Source,
        ))
        .id();

    let first_species = commands
        .spawn(Species {
            bodymass: settings.mean_bodymass_ratio_predator_prey,
            feeding_center: 0.0,
            feeding_range: sampler
                .0
                .gen_range(settings.feeding_range.min..=settings.feeding_range.max),
            first_occurence: 0.0,
            predator_strength: 1.0,
            dispersal_rate: 1,
        })
        .id();

    for x in 0..settings.grid_size.x {
        for y in 0..settings.grid_size.y {
            commands.spawn(Foodweb).with_children(|foodweb| {
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
    individuums: Query<&Individual>,
    species: Query<&Species>,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    let mut total_dispersal_rate = 0;
    for individuum in &individuums {
        total_dispersal_rate += species.get(individuum.0).unwrap().dispersal_rate;
    }

    let total_speciation_rate =
        settings.speciation_rate_per_individuum * individuums.iter().count() as u64;

    let increment = (1.0 + total_dispersal_rate as f64 / total_speciation_rate as f64)
        / settings.grid_size.x as f64
        / settings.grid_size.y as f64;
    time.0 += 1.0 / increment;

    let event_is_speciation = (total_dispersal_rate + total_speciation_rate)
        * sampler.0.gen_range(0..=1)
        >= total_dispersal_rate;

    if event_is_speciation {
        next_state.set(SimulationState::Speciation);
    } else {
        next_state.set(SimulationState::Dispersal);
    }
}

fn handle_speciation(
    time: Res<SimulationTime>,
    settings: Res<SimulationSettings>,
    mut sampler: ResMut<Sampler>,
    individuals: Query<(&Individual, &Parent), Without<Source>>,
    species: Query<&Species>,
    mut app_exit_events: ResMut<Events<AppExit>>,
    mut events: EventWriter<EventType>,
    mut commands: Commands,
) {
    let random = individuals.iter().choose(&mut sampler.0);
    if let Some((individual, parent)) = random {
        let og_species = species.get(individual.0).unwrap();
        let new = og_species.speciate(time.0, &settings, &mut sampler.0);
        if new.is_valid() {
            let new_id = commands.spawn(new).id();
            events.send(EventType::Speciation(parent.get(), new_id));
        }
    } else {
        info!("All species died");
        app_exit_events.send(AppExit::Success);
    }
}

fn handle_new_species(mut events: EventReader<EventType>) {
    for event in events.read() {
        if let EventType::Speciation(foodweb, species) = event {

        }
    }
}

fn handle_dispersal(mut next_state: ResMut<NextState<SimulationState>>) {
    next_state.set(SimulationState::Main);
}
