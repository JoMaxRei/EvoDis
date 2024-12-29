use std::collections::HashMap;

use bevy::prelude::*;
use rand::seq::SliceRandom;

use crate::*;

pub fn plugin(app: &mut App) {
    app.add_event::<SpeciationEvent>()
        .register_type::<FeedsOn>()
        .register_type::<PredationPressure>()
        .register_type::<Fitness>()
        .add_systems(OnEnter(SimulationState::Speciation), on_enter)
        .add_systems(
            Update,
            (
                handle_speciation,
                check_for_prey,
                // calculate,
                // calculate_predation_pressure,
                // calculate_gains,
                // check_survivability,
                // kill_species_with_no_individuals,
                next_state,
            )
                .chain()
                .run_if(in_state(SimulationState::Speciation)),
        )
        .add_systems(
            Update,
            (
                calculate,
                calculate_predation_pressure,
                calculate_gains,
                check_survivability,
            )
                .chain()
                .run_if(in_state(SimulationState::CheckSurvivors)),
        )
        .add_systems(
            OnExit(SimulationState::CheckSurvivors),
            (kill_species_with_no_individuals, on_exit),
        );
}

/// Contains all prey this entity feeds on
///
/// AKA num_preys and epsilon
#[derive(Component, Reflect)]
#[reflect(Component)]
struct FeedsOn(Vec<(Entity, f64)>);

#[derive(Event)]
enum SpeciationEvent {
    Insert(Entity),
}

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

#[derive(Resource, Reflect)]
#[reflect(Resource)]
struct NewestIndividual(Entity);

fn on_enter() {
    info!("handling speciations");
}

fn handle_speciation(
    time: Res<SimulationTime>,
    settings: Res<SimulationSettings>,
    mut sampler: ResMut<Sampler>,
    individuals: Query<(&Individual, &Parent), Without<Source>>,
    species: Query<&Species>,
    mut commands: Commands,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
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
        }
    } else {
        info!("All species died");
        next_state.set(SimulationState::Done);
    }
}

fn check_for_prey(
    speciation: Res<Speciation>,
    mut commands: Commands,
    individuals: Query<(&Individual, &Parent)>,
    species: Query<&Species>,
) {
            let new_species = species.get(speciation.species).unwrap();
            if individuals
                .iter()
                .filter(|(_, parent)| parent.get() == speciation.foodweb)
                .map(|(individiual, _)| species.get(individiual.0).unwrap())
                .any(|species| new_species.can_feed_on(species))
            {
                let newest_individual = commands
                    .spawn(Individual(speciation.species))
                    .set_parent(speciation.foodweb)
                    .id();
                commands.insert_resource(NewestIndividual(newest_individual));
                debug!(
                    "Added species {} as individual {} to foodweb {}",
                    speciation.species, newest_individual, speciation.foodweb
                );
            } else {
                debug!(
                    "Species {} has no prey in foodweb {}",
                    speciation.species, speciation.foodweb
                );
            }
}

fn calculate(
    speciation: Res<Speciation>,
    individuals: Query<(Entity, &Individual, &Parent)>,
    species: Query<&Species>,
    mut commands: Commands,
) {
            let mut inhabitants: Vec<(Entity, &Species)> = individuals
                .iter()
                .filter(|(_, _, parent)| parent.get() == speciation.foodweb)
                .map(|(entity, individuals, _)| (entity, species.get(individuals.0).unwrap()))
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
        
}

fn calculate_predation_pressure(
    predators: Query<&FeedsOn>,
    mut preys: Query<&mut PredationPressure>,
) {
    debug!("Start calculating predation pressure");
    for predator in &predators {
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
    inhabitants: Query<(Entity, &FeedsOn, Option<&Source>, &Individual)>,
    preys: Query<&PredationPressure>,
    species: Query<&Species>,
    mut commands: Commands,
) {
    debug!("Start calculating fitness");
    let mut sorted = inhabitants
        .iter()
        .map(|(e, f, o, i)| (e, f, o, species.get(i.0).unwrap()))
        .collect::<Vec<(Entity, &FeedsOn, Option<&Source>, &Species)>>();
    sorted.sort_by(|(_, _, _, lhs), (_, _, _, rhs)| lhs.bodymass.total_cmp(&rhs.bodymass));
    // E_Tilde
    let mut gains_without_competition = HashMap::new();
    // E
    let mut gains_with_competition = HashMap::new();
    for (entity, feeds_on, source_option, _) in sorted {
        if let Some(_) = source_option {
            gains_without_competition.insert(entity, settings.base_gain);
            gains_with_competition.insert(entity, settings.base_gain);
            continue;
        }

        for prey in &feeds_on.0 {
            let Some(&value) = gains_without_competition.get(&prey.0) else {
                warn!("There is no e_tilde for {}", prey.0);
                continue;
            };
            let predation_pressure = preys.get(prey.0).unwrap();
            debug!(
                "for {} e_tilde is {}, predation_loss is {}, epsilon is {}, sigma is {}",
                entity,
                value,
                predation_pressure.predation_loss(&settings),
                prey.1,
                predation_pressure.0
            );
            let gain = value * predation_pressure.predation_loss(&settings) * prey.1
                / predation_pressure.0;
            debug!("gain add for {} is {} from {}", entity, gain, prey.0);
            gains_without_competition
                .entry(entity)
                .and_modify(|e| *e += gain)
                .or_insert(gain);
            let value =
                gain / (1.0 + settings.competition_parameter * (predation_pressure.0 - prey.1));
            debug!("e add for {} is {} from {}", entity, value, prey.0);
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

    for (entity, _, _, _) in &inhabitants {
        let predation_pressure = preys.get(entity).unwrap();
        let value = 1.0
            / (predation_pressure.predation_loss(&settings)
                + 1.0 / gains_with_competition[&entity]);
        debug!(
            "for {} predation loss is {} and gains with comp is {}",
            entity,
            predation_pressure.predation_loss(&settings),
            gains_with_competition[&entity]
        );
        commands.entity(entity).insert(Fitness(value));
    }
    debug!("Finished calculating fitness");
}

fn check_survivability(
    newest_individual: Option<Res<NewestIndividual>>,
    mut sampler: ResMut<Sampler>,
    query: Query<(Entity, &Fitness), Without<Source>>,
    mut commands: Commands,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    commands.remove_resource::<NewestIndividual>();
    if let Some(newest) = newest_individual {
        let Ok((ent, new)) = query.get(newest.0) else {
            next_state.set(SimulationState::Done);
            return;
        };
        debug!("We have a newest {} with fitness {}", ent, new.0);
        if new.0 < 1.0 {
            commands.entity(ent).despawn_recursive();
            debug!("Killing newest");
            next_state.set(SimulationState::Main);
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
        debug!("Killing {}", dying);
        commands.entity(*dying).despawn_recursive();
        return;
    }
    next_state.set(SimulationState::Main);
}

fn kill_species_with_no_individuals(
    species_query: Query<Entity, (With<Species>, Without<Source>)>,
    individuals: Query<&Individual, Without<Source>>,
    mut commands: Commands,
) {
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
        }
    }
}

fn next_state(
    mut next_state: ResMut<NextState<SimulationState>>,
) {
        next_state.set(SimulationState::CheckSurvivors);
}

fn on_exit() {
    info!("done handling speciation");
}
