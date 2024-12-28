use bevy::prelude::*;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, evo_dis::plugin))
        .run();
}
