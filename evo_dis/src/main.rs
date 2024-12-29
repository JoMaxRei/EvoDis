use bevy::{
    log::{Level, LogPlugin}, prelude::*, state::app::StatesPlugin, window::PresentMode
};
use bevy_inspector_egui::quick::WorldInspectorPlugin;

fn main() {
    App::new()
        .add_plugins((
            MinimalPlugins,
            LogPlugin::default(),
            StatesPlugin::default(),
            // DefaultPlugins
            //     .set(LogPlugin {
            //         level: Level::INFO,
            //         ..default()
            //     })
            //     .set(WindowPlugin {
            //         primary_window: Some(Window {
            //             present_mode: PresentMode::AutoNoVsync,
            //             ..default()
            //         }),
            //         ..default()
            //     }),
            // WorldInspectorPlugin::new(),
            evo_dis::plugin,
        ))
        .run();
}
