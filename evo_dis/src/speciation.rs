use bevy::prelude::*;

use crate::*;

pub fn plugin(app: &mut App) {
    app.add_systems(
        Update,
        (handle_speciation, next_state)
            .chain()
            .run_if(in_state(SimulationState::Speciation)),
    );
}

fn handle_speciation(
    mut statistics: ResMut<Statistics>,
    time: Res<SimulationTime>,
    settings: Res<SimulationSettings>,
    mut sampler: ResMut<Sampler>,
    individuals: Query<(&Individual, &Parent), Without<Source>>,
    species: Query<&Species>,
    mut commands: Commands,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    let _ = info_span!("handle_speciation", name = "handle_speciation").entered();
    let random = individuals.iter().choose(&mut sampler.0);
    if let Some((individual, parent)) = random {
        let og_species = species.get(individual.0).unwrap();
        let new = og_species.speciate(time.0, &settings, &mut sampler.0);
        if new.is_valid() {
            let new_id = commands
                .spawn((Name::from(format!("Species {}", time.0)), new))
                .id();
            commands.insert_resource(Speciation {
                foodweb: parent.get(),
                species: new_id,
            });
            statistics.speciations += 1;
        }
    } else {
        info!("All species died");
        next_state.set(SimulationState::Done);
    }
}

fn next_state(mut next_state: ResMut<NextState<SimulationState>>) {
    next_state.set(SimulationState::Prey);
}
