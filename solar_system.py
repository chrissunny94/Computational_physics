import pygame
import math
import sys

pygame.init()

# ---------------- CONFIG ---------------- #
WIDTH, HEIGHT = 1200, 800
WIN = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Solar System Simulation")

FONT = pygame.font.SysFont("Arial", 16)

WHITE = (255, 255, 255)
YELLOW = (255, 255, 0)
BLUE = (100, 149, 237)
RED = (188, 39, 50)
GREY = (130, 130, 130)

# Scale: 1 AU = 150 pixels (adjust for zoom)
SCALE = 150  
TIMESTEP = 60 * 60 * 24  # 1 day per frame


# -----------------------------------------
# PLANET CLASS
# -----------------------------------------
class Planet:
    AU = 1.496e11        # meters
    G = 6.674e-11        # gravitational constant

    def __init__(self, x, y, radius, color, mass):
        self.x = x            # actual position in meters
        self.y = y
        self.radius = radius  # for drawing
        self.color = color
        self.mass = mass

        self.orbit = []
        self.sun = False
        self.distance_to_sun = 0

        self.x_vel = 0
        self.y_vel = 0

    # -------------------------
    # GRAVITY + MOTION
    # -------------------------
    def attraction(self, other):
        dx = other.x - self.x
        dy = other.y - self.y
        distance = math.sqrt(dx ** 2 + dy ** 2)

        if other.sun:
            self.distance_to_sun = distance

        force = self.G * self.mass * other.mass / distance ** 2
        angle = math.atan2(dy, dx)
        fx = math.cos(angle) * force
        fy = math.sin(angle) * force
        return fx, fy

    def update_position(self, planets):
        total_fx = total_fy = 0
        for planet in planets:
            if planet == self:
                continue
            fx, fy = self.attraction(planet)
            total_fx += fx
            total_fy += fy

        # Update velocities
        self.x_vel += total_fx / self.mass * TIMESTEP
        self.y_vel += total_fy / self.mass * TIMESTEP

        # Update positions
        self.x += self.x_vel * TIMESTEP
        self.y += self.y_vel * TIMESTEP

        # Save orbit path
        self.orbit.append((self.x, self.y))
        if len(self.orbit) > 400:
            self.orbit.pop(0)

    # -------------------------
    # DRAWING
    # -------------------------
    def draw(self, win):
        # Convert from meters → screen pixels
        x = self.x / self.AU * SCALE + WIDTH / 2
        y = self.y / self.AU * SCALE + HEIGHT / 2

        # Draw orbit path
        if len(self.orbit) > 2:
            points = []
            for px, py in self.orbit:
                sx = px / self.AU * SCALE + WIDTH / 2
                sy = py / self.AU * SCALE + HEIGHT / 2
                points.append((sx, sy))
            pygame.draw.lines(win, self.color, False, points, 2)

        # Draw planet body
        pygame.draw.circle(win, self.color, (int(x), int(y)), self.radius)

        # Show distance label
        if not self.sun:
            dist_text = FONT.render(f"{self.distance_to_sun/1e9:.1f} Gm", True, WHITE)
            win.blit(dist_text, (int(x) + 10, int(y) - 10))


# ---------------------------
# CREATE PLANETS
# ---------------------------
def make_planets():
    sun = Planet(0, 0, 30, YELLOW, 1.98892e30)
    sun.sun = True

    earth = Planet(-Planet.AU, 0, 10, BLUE, 5.9742e24)
    earth.y_vel = 29.783e3  # Earth's orbital velocity

    mars = Planet(-1.524 * Planet.AU, 0, 8, RED, 6.39e23)
    mars.y_vel = 24.077e3

    mercury = Planet(0.387 * Planet.AU, 0, 6, GREY, 3.30e23)
    mercury.y_vel = -47.4e3

    venus = Planet(0.723 * Planet.AU, 0, 8, (255, 200, 0), 4.8685e24)
    venus.y_vel = -35.02e3

    return [sun, earth, mars, mercury, venus]


# ---------------------------
# MAIN LOOP
# ---------------------------
def main():
    clock = pygame.time.Clock()
    running = True
    planets = make_planets()

    while running:
        clock.tick(60)
        WIN.fill((0, 0, 0))

        # Quit
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        # Update + draw planets
        for planet in planets:
            planet.update_position(planets)
            planet.draw(WIN)

        pygame.display.update()

    pygame.quit()
    sys.exit()


if __name__ == "__main__":
    main()
