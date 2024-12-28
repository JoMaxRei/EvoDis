use bevy::prelude::*;
use rand::Rng;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Distribution, StandardNormal};

use crate::settings::SimulationSettings;

#[derive(Component)]
pub struct Species {
    pub bodymass: f64,
    pub feeding_center: f64,
    pub feeding_range: f64,
    pub first_occurence: f64,
    pub predator_strength: f64,
    pub dispersal_rate: u64,
}

impl Species {
    pub fn speciate(
        &self,
        first_occurence: f64,
        settings: &SimulationSettings,
        sampler: &mut ChaCha8Rng,
    ) -> Self {
        let bodymass = self.bodymass
            + sampler.gen_range(settings.bodymass_range.min..=settings.bodymass_range.max);
        let mut feeding_center = 0.0;
        let normal = StandardNormal;
        let mut feeding_range = 0.0;
        loop {
            loop {
                let value: f64 = normal.sample(sampler);
                feeding_center = bodymass - settings.mean_bodymass_ratio_predator_prey + value;
                if feeding_center - settings.feeding_range.max <= bodymass {
                    break;
                }
            }

            feeding_range =
                sampler.gen_range(settings.feeding_range.min..=settings.feeding_range.max);

            if feeding_center - feeding_range <= bodymass {
                break;
            }
        }
        let min = (self.dispersal_rate as f64 / (1.0 + settings.dispersel_variance)).log2();
        let diff = (self.dispersal_rate as f64 * (1.0 + settings.dispersel_variance)).log2() - min;
        let dispersal_rate = 2_f64.powf(sampler.r#gen::<f64>() * diff + min).round() as u64;

        Species {
            bodymass,
            feeding_center,
            feeding_range,
            first_occurence,
            predator_strength: Self::calculate_predator_strength(
                dispersal_rate,
                settings.zero_crossing,
                settings.initial_dispersal_rate,
            ),
            dispersal_rate,
        }
    }

    fn calculate_predator_strength(
        dispersal_rate: u64,
        zero_crossing: f64,
        initial_dispersal_rate: f64,
    ) -> f64 {
        ((dispersal_rate as f64 / zero_crossing).log10()
            / (initial_dispersal_rate / zero_crossing).log10())
        .clamp(0.0, 1.0)
    }

    pub fn is_valid(&self) -> bool {
        if self.bodymass <= 0.0 {
            return false;
        }

        true
    }
}

#[derive(Component)]
pub struct Individual(pub Entity);

#[derive(Component)]
pub struct Source;
