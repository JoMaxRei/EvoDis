use std::collections::HashMap;

use bevy::prelude::*;
use rand::seq::SliceRandom;

use crate::{
    Sampler, SimulationState, Speciation,
    foodweb::Foodweb,
    settings::{SimulationSettings, Statistics},
    species::{Individual, Source, Species},
};

pub fn plugin(app: &mut App) {
    app.register_type::<FeedsOn>()
        .register_type::<PredationPressure>()
        .register_type::<Fitness>()
        .add_systems(OnEnter(SimulationState::Prey), tag_individuals)
        .add_systems(
            Update,
            check_for_prey.run_if(in_state(SimulationState::Prey)),
        )
        .add_systems(OnEnter(SimulationState::Main), remove_tags)
        .add_systems(
            Update,
            (
                calculate_foodweb,
                calculate_predation_pressure,
                calculate_gains,
                check_survivability,
            )
                .chain()
                .run_if(in_state(SimulationState::CheckSurvivors)),
        )
        .add_systems(
            OnExit(SimulationState::CheckSurvivors),
            kill_species_with_no_individuals,
        );
}

#[derive(Resource, Reflect)]
#[reflect(Resource)]
struct NewestIndividual(Entity);

/// Contains all prey this entity feeds on
///
/// AKA num_preys and epsilon
#[derive(Component, Reflect)]
#[reflect(Component)]
struct FeedsOn(Vec<(Entity, f64)>);

/// Contains sigma
#[derive(Component, Reflect)]
#[reflect(Component)]
struct PredationPressure(f64);

#[derive(Component, Reflect)]
#[reflect(Component)]
struct Fitness(f64);

impl PredationPressure {
    /// AKA P
    fn predation_loss(&self, settings: &SimulationSettings) -> f64 {
        settings.predation_parameter * self.0 / (settings.predation_parameter * self.0 + 1.0)
    }
}

#[derive(Component)]
struct CurrentIndividual;

fn tag_individuals(
    speciation: Res<Speciation>,
    foodwebs: Query<&Children, With<Foodweb>>,
    mut commands: Commands,
) {
    let _ = info_span!("tag_individuals", name = "tag_individuals").entered();
    if let Ok(foodweb) = foodwebs.get(speciation.foodweb) {
        for child in foodweb {
            commands.entity(*child).insert(CurrentIndividual);
        }
        debug!("tagged {} children before adding new", foodweb.len());
    }
}

fn check_for_prey(
    speciation: Res<Speciation>,
    mut commands: Commands,
    individuals: Query<&Individual, With<CurrentIndividual>>,
    species: Query<&Species>,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    let _ = info_span!("check_for_prey", name = "check_for_prey").entered();
    debug!("Start checking for prey");
    let new_species = species.get(speciation.species).unwrap();
    if individuals
        .iter()
        .map(|individiual| species.get(individiual.0).unwrap())
        .any(|species| new_species.can_feed_on(species))
    {
        let newest_individual = commands
            .spawn((Individual(speciation.species), CurrentIndividual))
            .set_parent(speciation.foodweb)
            .id();
        commands.insert_resource(NewestIndividual(newest_individual));
        trace!(
            "Added species {} as individual {} to foodweb {}",
            speciation.species, newest_individual, speciation.foodweb
        );
        next_state.set(SimulationState::CheckSurvivors);
    } else {
        trace!(
            "Species {} has no prey in foodweb {}",
            speciation.species, speciation.foodweb
        );
        next_state.set(SimulationState::Main);
    }
    debug!("Finished checking for prey");
}

fn calculate_foodweb(
    individuals: Query<(Entity, &Individual), With<CurrentIndividual>>,
    species: Query<&Species>,
    mut commands: Commands,
) {
    let _ = info_span!("calculate_foodweb", name = "calculate_foodweb").entered();
    debug!("Start calculating foodweb");
    let mut inhabitants: Vec<(Entity, &Species)> = individuals
        .iter()
        .map(|(entity, individuals)| (entity, species.get(individuals.0).unwrap()))
        .collect();
    inhabitants.sort_by(|(_, lhs), (_, rhs)| lhs.bodymass.total_cmp(&rhs.bodymass));

    for index in 0..inhabitants.len() {
        let (entity, species) = inhabitants[index];
        let mut prey = Vec::new();
        for compare_index in 0..index {
            let (prey_entity, rhs) = inhabitants[compare_index];
            if species.can_feed_on(rhs) {
                prey.push((prey_entity, species.epsilon(rhs)));
            }
        }
        commands
            .entity(entity)
            .insert((FeedsOn(prey), PredationPressure(0.0)));
    }
    debug!("Finished calculating foodweb");
}

fn calculate_predation_pressure(
    speciation: Res<Speciation>,
    predators: Query<(Entity, &FeedsOn), With<CurrentIndividual>>,
    mut preys: Query<&mut PredationPressure>,
) {
    let _ = info_span!(
        "calculate_predation_pressure",
        name = "calculate_predation_pressure"
    )
    .entered();
    debug!("Start calculating predation pressure");
    for (entity, predator) in &predators {
        trace!(
            "predator {:?} on {} has {} preys",
            entity,
            speciation.foodweb,
            predator.0.len()
        );
        for (entity, pressure) in &predator.0 {
            if let Ok(mut prey) = preys.get_mut(*entity) {
                prey.0 += pressure;
            }
        }
    }
    debug!("Finished calculating predation pressure");
}

fn calculate_gains(
    settings: Res<SimulationSettings>,
    individuals: Query<(Entity, &FeedsOn, Option<&Source>, &Individual), With<CurrentIndividual>>,
    preys: Query<&PredationPressure>,
    species: Query<&Species>,
    mut commands: Commands,
) {
    let _ = info_span!("calculate_gains", name = "calculate_gains").entered();
    debug!("Start calculating fitness");
    let mut inhabitants = individuals
        .iter()
        .map(|(e, f, o, i)| (e, f, o, species.get(i.0).unwrap()))
        .collect::<Vec<(Entity, &FeedsOn, Option<&Source>, &Species)>>();
    inhabitants.sort_by(|(_, _, _, lhs), (_, _, _, rhs)| lhs.bodymass.total_cmp(&rhs.bodymass));
    // E_Tilde
    let mut gains_without_competition = HashMap::new();
    // E
    let mut gains_with_competition = HashMap::new();
    for (entity, feeds_on, source_option, _) in inhabitants {
        if let Some(_) = source_option {
            gains_without_competition.insert(entity, settings.base_gain);
            gains_with_competition.insert(entity, settings.base_gain);
            continue;
        }

        trace!("{} can feed on {:?}", entity, feeds_on.0);
        for prey in &feeds_on.0 {
            let Some(&value) = gains_without_competition.get(&prey.0) else {
                warn!("There is no e_tilde for {}", prey.0);
                continue;
            };
            let predation_pressure = preys.get(prey.0).unwrap();
            trace!(
                "for {} e_tilde is {}, predation_loss is {}, epsilon is {}, sigma is {}",
                entity,
                value,
                predation_pressure.predation_loss(&settings),
                prey.1,
                predation_pressure.0
            );
            let gain = value * predation_pressure.predation_loss(&settings) * prey.1
                / predation_pressure.0;
            trace!("gain add for {} is {} from {}", entity, gain, prey.0);
            gains_without_competition
                .entry(entity)
                .and_modify(|e| *e += gain)
                .or_insert(gain);
            let value =
                gain / (1.0 + settings.competition_parameter * (predation_pressure.0 - prey.1));
            trace!("e add for {} is {} from {}", entity, value, prey.0);
            gains_with_competition
                .entry(entity)
                .and_modify(|e| *e += value)
                .or_insert(value);
        }

        gains_without_competition
            .entry(entity)
            .and_modify(|e| *e *= settings.trophic_loss);
        gains_with_competition
            .entry(entity)
            .and_modify(|e| *e *= settings.trophic_loss);
    }

    for (entity, _, _, _) in &individuals {
        let predation_pressure = preys.get(entity).unwrap();
        let gains = gains_with_competition.get(&entity).unwrap_or(&0.0);
        let value = 1.0 / (predation_pressure.predation_loss(&settings) + 1.0 / gains);
        trace!(
            "for {} predation loss is {} and gains with comp is {}, total fitness is {}",
            entity,
            predation_pressure.predation_loss(&settings),
            gains,
            value,
        );
        commands.entity(entity).insert(Fitness(value));
    }
    debug!("Finished calculating fitness");
}

fn check_survivability(
    mut statistics: ResMut<Statistics>,
    newest_individual: Option<Res<NewestIndividual>>,
    mut sampler: ResMut<Sampler>,
    query: Query<(Entity, &Fitness), (Without<Source>, With<CurrentIndividual>)>,
    mut commands: Commands,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    let _ = info_span!("check_survivability", name = "check_survivability").entered();
    debug!("Start calculating survivability");
    commands.remove_resource::<NewestIndividual>();
    if let Some(newest) = newest_individual {
        let Ok((ent, new)) = query.get(newest.0) else {
            next_state.set(SimulationState::Wait);
            return;
        };
        trace!("We have a newest {} with fitness {}", ent, new.0);
        if new.0 < 1.0 {
            commands.entity(ent).despawn_recursive();
            next_state.set(SimulationState::Main);
            debug!(
                "Finished calculating survivability with killing newest {}",
                ent
            );
            statistics.died_individuals += 1;
            return;
        }
    }
    let mut fitnesses = query
        .iter()
        .filter(|(_, fitness)| fitness.0 < 1.0)
        .collect::<Vec<(Entity, &Fitness)>>();
    fitnesses.sort_by(|(_, lhs), (_, rhs)| lhs.0.total_cmp(&rhs.0));
    if let Some((_, first_fitness)) = fitnesses.first() {
        let contenders = fitnesses
            .iter()
            .filter(|(_, element)| element.0 <= first_fitness.0)
            .collect::<Vec<&(Entity, &Fitness)>>();
        let (dying, _) = contenders.choose(&mut sampler.0).unwrap();
        commands.entity(*dying).despawn_recursive();
        debug!("Finished calculating survivability with killing {}", dying);
        statistics.died_individuals += 1;
        return;
    }
    next_state.set(SimulationState::Wait);
    debug!("Finished calculating survivability");
}

fn kill_species_with_no_individuals(
    mut statistics: ResMut<Statistics>,
    species_query: Query<Entity, (With<Species>, Without<Source>)>,
    individuals: Query<&Individual, Without<Source>>,
    mut commands: Commands,
) {
    let _ = info_span!(
        "kill_species_with_no_individuals",
        name = "kill_species_with_no_individuals"
    )
    .entered();
    for species in &species_query {
        if individuals
            .iter()
            .filter(|individual| individual.0 == species)
            .count()
            == 0
        {
            commands.entity(species).despawn_recursive();
            debug!(
                "Species {} does not have any individuals left, killing...",
                species
            );
            statistics.died_species += 1;
        }
    }
}

fn remove_tags(query: Query<Entity, With<CurrentIndividual>>, mut commands: Commands) {
    let _ = info_span!("remove_tags", name = "remove_tags").entered();
    for entity in &query {
        commands.entity(entity).remove::<CurrentIndividual>();
    }
    debug!("de-tagged {} children", query.iter().len());
}
