use bevy::{
    log::{Level, LogPlugin},
    prelude::*,
};
use bevy_inspector_egui::quick::WorldInspectorPlugin;

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins.set(LogPlugin {
                level: Level::DEBUG,
                ..default()
            }),
            WorldInspectorPlugin::new(),
            evo_dis::plugin,
        ))
        .run();
}
