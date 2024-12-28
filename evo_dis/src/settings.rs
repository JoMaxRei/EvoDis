use bevy::prelude::*;

pub struct Range {
    pub min: f64,
    pub max: f64,
}

impl Range {
    pub const fn new(min: f64, max: f64) -> Self {
        Range { min, max }
    }
}

#[derive(Resource)]
pub struct SimulationSettings {
    pub speciation_rate_per_individuum: u64,
    pub feeding_range: Range,
    pub grid_size: UVec2,
    pub mean_bodymass_ratio_predator_prey: f64,
    pub bodymass_range: Range,
    pub dispersel_variance: f64,
    pub zero_crossing: f64,
    pub initial_dispersal_rate: f64,
}

impl Default for SimulationSettings {
    fn default() -> Self {
        SimulationSettings {
            speciation_rate_per_individuum: 1_000_000,
            feeding_range: Range::new(0.2, 2.0),
            grid_size: UVec2::new(30, 30),
            mean_bodymass_ratio_predator_prey: 2.0,
            bodymass_range: Range::new(-5.0_f64.log10(), 5.0_f64.log10()),
            dispersel_variance: 0.3,
            zero_crossing: 10_000.0,
            initial_dispersal_rate: 10.0,
        }
    }
}
