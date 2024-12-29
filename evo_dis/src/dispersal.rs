use bevy::prelude::*;
use rand::Rng;
use rand_distr::{Distribution, WeightedIndex};

use crate::{
    Sampler, SimulationState, Speciation,
    foodweb::Foodweb,
    position::Position,
    settings::{SimulationSettings, Statistics},
    species::{Individual, Source, Species},
};

pub fn plugin(app: &mut App) {
    app.add_systems(
        Update,
        (find_dispersing_individual)
            .chain()
            .run_if(in_state(SimulationState::Dispersal)),
    );
}

fn find_dispersing_individual(
    mut statistics: ResMut<Statistics>,
    settings: Res<SimulationSettings>,
    mut sampler: ResMut<Sampler>,
    species: Query<&Species>,
    individuals: Query<(Entity, &Individual, &Parent), Without<Source>>,
    foodwebs: Query<(Entity, &Position), With<Foodweb>>,
    mut next_state: ResMut<NextState<SimulationState>>,
    mut commands: Commands,
) {
    statistics.dispersals += 1;
    let all = individuals
        .iter()
        .map(|(entity, individual, parent)| {
            (
                entity,
                species.get(individual.0).unwrap().dispersal_rate,
                parent,
            )
        })
        .collect::<Vec<(Entity, u64, &Parent)>>();
    let weights = all.iter().map(|(_, disp, _)| *disp).collect::<Vec<u64>>();
    let dist = WeightedIndex::new(&weights).unwrap();
    let index = dist.sample(&mut sampler.0);
    let (entity, _, parent) = all[index];
    let (_, position) = foodwebs.get(parent.get()).unwrap();

    let (x, y) = find_neigbhor(sampler.0.gen_range(0.0..1.0), &settings, position);

    let (_, original_individual, _) = individuals.get(entity).unwrap();
    if let Some((foodweb, _)) = foodwebs
        .iter()
        .find(|(_, web)| **web == Position::new(x, y))
    {
        if individuals.iter().any(|(_, individual, parent)| {
            parent.get() == foodweb && individual == original_individual
        }) {
            trace!(
                "Habitat {} already has the species {} coming from {}",
                foodweb,
                original_individual.0,
                parent.get()
            );
            next_state.set(SimulationState::Main);
            statistics.failed_dispersals += 1;
        } else if individuals
            .iter()
            .any(|(_, _, parent)| parent.get() == foodweb)
        {
            commands.insert_resource(Speciation {
                foodweb,
                species: original_individual.0,
            });
            next_state.set(SimulationState::Prey);
        }
    }
}

fn find_neigbhor(
    random_value: f64,
    settings: &SimulationSettings,
    position: &Position,
) -> (u32, u32) {
    let kernel_size = 2 * settings.range + 1;
    let kernel_size_sqared = kernel_size.pow(2);
    let random = random_value * (kernel_size_sqared as f64 - 1.0);
    let source = ((kernel_size_sqared + 1) / 2 + random as u32) % kernel_size_sqared;
    let x = (settings.range * settings.grid_size.x + position.x + (source / kernel_size)
        - settings.range)
        % settings.grid_size.x;
    let y = (settings.range * settings.grid_size.y + position.x + (source % kernel_size)
        - settings.range)
        % settings.grid_size.y;
    (x, y)
}

#[test]
fn test_find_neighbor() {
    assert_eq!(
        find_neigbhor(
            0.5,
            &SimulationSettings {
                range: 1,
                grid_size: UVec2::new(1, 1),
                ..default()
            },
            &Position::new(0, 0)
        ),
        (0, 0)
    );
}
