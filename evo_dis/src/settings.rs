use bevy::prelude::*;

#[derive(Reflect)]
pub struct Range {
    pub min: f64,
    pub max: f64,
}

impl Range {
    pub const fn new(min: f64, max: f64) -> Self {
        Range { min, max }
    }
}

#[derive(Resource, Reflect)]
#[reflect(Resource)]
pub struct SimulationSettings {
    pub speciation_rate_per_individuum: u64,
    pub feeding_range: Range,
    pub grid_size: UVec2,
    pub mean_bodymass_ratio_predator_prey: f64,
    pub bodymass_range: Range,
    pub dispersel_variance: f64,
    pub zero_crossing: f64,
    pub initial_dispersal_rate: f64,
    /// AKA p
    pub predation_parameter: f64,
    /// AKA E0
    pub base_gain: f64,
    /// AKA xi
    pub competition_parameter: f64,
    /// AKA x
    pub trophic_loss: f64,
}

impl Default for SimulationSettings {
    fn default() -> Self {
        let base_rate = 1_000_000;
        SimulationSettings {
            speciation_rate_per_individuum: base_rate,
            feeding_range: Range::new(0.2, 2.0),
            grid_size: UVec2::new(1, 1),
            mean_bodymass_ratio_predator_prey: 2.0,
            bodymass_range: Range::new(-5.0_f64.log10(), 5.0_f64.log10()),
            dispersel_variance: 0.3,
            zero_crossing: 10_000.0 * base_rate as f64,
            initial_dispersal_rate: 10.0 * base_rate as f64,
            predation_parameter: 1.0,
            base_gain: 1000.0,
            competition_parameter: 3.0,
            trophic_loss: 0.4,
        }
    }
}
