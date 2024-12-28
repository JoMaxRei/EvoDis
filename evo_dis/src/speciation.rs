use std::collections::HashMap;

use bevy::prelude::*;

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
                calculate,
                calculate_predation_pressure,
                calculate_gains,
                next_state,
            )
                .chain()
                .run_if(in_state(SimulationState::Speciation)),
        )
        .add_systems(OnExit(SimulationState::Speciation), on_exit);
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

fn on_enter() {
    info!("handling speciations");
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
            let new_id = commands
                .spawn((Name::from(format!("Species {}", time.0)), new))
                .id();
            events.send(EventType::Speciation(parent.get(), new_id));
        }
    } else {
        info!("All species died");
        app_exit_events.send(AppExit::Success);
    }
}

fn check_for_prey(
    mut events: EventReader<EventType>,
    mut commands: Commands,
    individuals: Query<(&Individual, &Parent)>,
    species: Query<&Species>,
    mut speciation_events: EventWriter<SpeciationEvent>,
) {
    for event in events.read() {
        if let EventType::Speciation(foodweb_entity, species_entity) = event {
            let new_species = species.get(*species_entity).unwrap();
            if individuals
                .iter()
                .filter(|(_, parent)| parent.get() == *foodweb_entity)
                .map(|(individiual, _)| species.get(individiual.0).unwrap())
                .any(|species| new_species.can_feed_on(species))
            {
                commands
                    .entity(*foodweb_entity)
                    .with_child(Individual(*species_entity));
                info!(
                    "Added species {} to foodweb {}",
                    species_entity, foodweb_entity
                );
                speciation_events.send(SpeciationEvent::Insert(*foodweb_entity));
            } else {
                warn!(
                    "Species {} is not viable in foodweb {}",
                    species_entity, foodweb_entity
                );
            }
        }
    }
}

fn calculate(
    mut events: EventReader<SpeciationEvent>,
    individuals: Query<(Entity, &Individual, &Parent)>,
    species: Query<&Species>,
    mut commands: Commands,
) {
    for event in events.read() {
        if let SpeciationEvent::Insert(foodweb_entity) = event {
            let mut inhabitants: Vec<(Entity, &Species)> = individuals
                .iter()
                .filter(|(_, _, parent)| parent.get() == *foodweb_entity)
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
                panic!("There is no e_tilde for {}", prey.0);
            };
            let predation_pressure = preys.get(prey.0).unwrap();
            info!("for {} e_tilde is {}, predation_loss is {}, epsilon is {}, sigma is {}", entity, value, predation_pressure.predation_loss(&settings), prey.1, predation_pressure.0);
            let gain = value * predation_pressure.predation_loss(&settings) * prey.1
                / predation_pressure.0;
            info!("gain add for {} is {} from {}", entity, gain, prey.0);
            gains_without_competition
                .entry(entity)
                .and_modify(|e| *e += gain)
                .or_insert(gain);
            let value =
                gain / (1.0 + settings.competition_parameter * (predation_pressure.0 - prey.1));
            info!("e add for {} is {} from {}", entity, value, prey.0);
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
        info!(
            "for {} predation loss is {} and gains with comp is {}",
            entity,
            predation_pressure.predation_loss(&settings),
            gains_with_competition[&entity]
        );
        commands.entity(entity).insert(Fitness(value));
    }
}

fn next_state(mut next_state: ResMut<NextState<SimulationState>>) {
    next_state.set(SimulationState::Done);
}

fn on_exit() {
    info!("done handling speciation");
}
