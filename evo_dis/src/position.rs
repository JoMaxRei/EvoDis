use bevy::prelude::*;

#[derive(Component, Reflect, PartialEq)]
#[reflect(Component)]
pub struct Position {
    pub x: u32,
    pub y: u32,
}

impl Position {
    pub const fn new(x: u32, y: u32) -> Self {
        Position { x, y }
    }
}
