#include <iostream>
#include <chrono>
#include <string>
#include <algorithm>
#include <math.h>
#include <vector>


int skillEffectsLENGTH = 0;
int unitsLENGTH = 0;

constexpr double MAP_RADIUS = 6000.0;
int TANKERS_BY_PLAYER;
constexpr int TANKERS_BY_PLAYER_MIN = 1;
constexpr int TANKERS_BY_PLAYER_MAX = 3;

constexpr double WATERTOWN_RADIUS = 3000.0;

constexpr int TANKER_THRUST = 500;
constexpr double TANKER_EMPTY_MASS = 2.5;
constexpr double TANKER_MASS_BY_WATER = 0.5;
constexpr double TANKER_FRICTION = 0.40;
constexpr double TANKER_RADIUS_BASE = 400.0;
constexpr double TANKER_RADIUS_BY_SIZE = 50.0;
constexpr int TANKER_EMPTY_WATER = 1;
constexpr int TANKER_MIN_SIZE = 4;
constexpr int TANKER_MAX_SIZE = 10;
constexpr double TANKER_MIN_RADIUS = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * TANKER_MIN_SIZE;
constexpr double TANKER_MAX_RADIUS = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * TANKER_MAX_SIZE;
constexpr double TANKER_SPAWN_RADIUS = 8000.0;
constexpr int TANKER_START_THRUST = 2000;

constexpr int MAX_THRUST = 300;
constexpr int MAX_RAGE = 300;
constexpr int WIN_SCORE = 50;

constexpr double REAPER_MASS = 0.5;
constexpr double REAPER_FRICTION = 0.20;
constexpr int REAPER_SKILL_DURATION = 3;
constexpr int REAPER_SKILL_COST = 30;
constexpr int REAPER_SKILL_ORDER = 0;
constexpr double REAPER_SKILL_RANGE = 2000.0;
constexpr double REAPER_SKILL_RADIUS = 1000.0;
constexpr double REAPER_SKILL_MASS_BONUS = 10.0;

constexpr double DESTROYER_MASS = 1.5;
constexpr double DESTROYER_FRICTION = 0.30;
constexpr int DESTROYER_SKILL_DURATION = 1;
constexpr int DESTROYER_SKILL_COST = 60;
constexpr int DESTROYER_SKILL_ORDER = 2;
constexpr double DESTROYER_SKILL_RANGE = 2000.0;
constexpr double DESTROYER_SKILL_RADIUS = 1000.0;
constexpr int DESTROYER_NITRO_GRENADE_POWER = 1000;

constexpr double DOOF_MASS = 1.0;
constexpr double DOOF_FRICTION = 0.25;
constexpr double DOOF_RAGE_COEF = 1.0 / 100.0;
constexpr int DOOF_SKILL_DURATION = 3;
constexpr int DOOF_SKILL_COST = 30;
constexpr int DOOF_SKILL_ORDER = 1;
constexpr double DOOF_SKILL_RANGE = 2000.0;
constexpr double DOOF_SKILL_RADIUS = 1000.0;

constexpr double LOOTER_RADIUS = 400.0;
constexpr int LOOTER_REAPER = 0;
constexpr int LOOTER_DESTROYER = 1;
constexpr int LOOTER_DOOF = 2;

constexpr int TYPE_TANKER = 3;
constexpr int TYPE_WRECK = 4;
constexpr int TYPE_REAPER_SKILL_EFFECT = 5;
constexpr int TYPE_DOOF_SKILL_EFFECT = 6;
constexpr int TYPE_DESTROYER_SKILL_EFFECT = 7;

constexpr double EPSILON = 0.00001;
constexpr double MIN_IMPULSE = 30.0;
constexpr double IMPULSE_COEFF = 0.5;

// Global first free id for all elements on the map 
int GLOBAL_ID = 0;

class Point;
class Wreck;
class Player;
class SkillEffect;
class Collision;
class Doof;
class Destroyer;
class Reaper;
class Unit;
class Looter;
class SkillEffect;
class ReaperSkillEffect;
class DoofSkillEffect;
class DestroyerSkillEffect;
class Tanker;

enum Skill {
    DOOF, DESTROYER, REAPER
};

enum Action {
    SKILL, MOVE, WAIT
};

// ****************************************************************************************

//TODO: might need to add equals and hash function for class Point

class Point {
public:
    double x;
    double y;

    Point() {};

    Point(double x, double y) {
        this->x = x;
        this->y = y;
    }
};

inline double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

inline double dist2(double x1, double y1, double x2, double y2) {
    return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
}

inline double dist(Point* p, double x2, double y2) {
    return dist(p->x, p->y, x2, y2);
}

inline double dist2(Point* p, double x2, double y2) {
    return dist2(p->x, p->y, x2, y2);
}

inline double dist(Point* u1, Point* u2) {
    return dist(u1->x, u1->y, u2->x, u2->y);
}

inline double dist2(Point* u1, Point* u2) {
    return dist2(u1->x, u1->y, u2->x, u2->y);
}

void moveTo(Point* u1, Point* u2, double distance) {
    double d = dist(u1, u2);

    if (d < EPSILON) {
        return;
    }

    double dx = u1->x - u2->x;
    double dy = u1->y - u2->y;
    double coef = distance / d;

    u1->x += dx * coef;
    u1->y += dy * coef;
}

inline bool isInRange(Point* p1, Point* p2, double range) {
    return p1 != p2 && dist(p1, p2) <= range;
}

// ****************************************************************************************
class Collision {
public:
    double t;
    Unit* a;
    Unit* b;

    Collision(double t) : t(t) , a(nullptr) , b(nullptr) {

    }

    Collision(double t, Unit* a) : t(t), a(a), b(nullptr) {

    }

    Collision(double t, Unit* a, Unit* b) {
        this->t = t;
        this->a = a;
        this->b = b;
    }

    Tanker* dead();
};

// ****************************************************************************************

// Center of the map
Point WATERTOWN(0, 0);

Collision NULL_COLLISION(1.0 + EPSILON);

// ****************************************************************************************
class Player {
public:
    int score;
    int index;
    int rage;
    Looter* looters[3];
    bool dead;

    Player(int index) {
        this->index = index;
    }

    void kill() {
        dead = true;
    }

    Reaper* getReaper();

    Destroyer* getDestroyer();

    Doof* getDoof();
};
// ****************************************************************************************
class Wreck : public Point {
public:
    int id;
    double radius;
    int water;
    bool known;
    Player* player;

    Wreck(double x, double y, int water, double radius) : Point(x, y) {
        id = GLOBAL_ID++;
        this->radius = radius;
        this->water = water;
    }

    bool harvest(Player* players[], std::vector<SkillEffect*>& skillEffects);

};

// ****************************************************************************************

//TODO: may need to write equals and hash functions
class Unit : public Point {
public:
    int type;
    int id;
    double vx;
    double vy;
    double radius;
    double mass;
    double friction;
    bool known;

    Unit(int type , double x, double y) : Point(x, y) {
        id = GLOBAL_ID++;
        this->type = type;
        vx = 0.0;
        vy = 0.0;
        known = false;
    }

    void move(double t) {
        x += vx * t;
        y += vy * t;
    }

    double speed() {
        return sqrt(vx * vx + vy * vy);
    }

    void thrust(Point* p, int power) {
        double distance = dist(this, p);

        if (distance <= EPSILON)
            return;

        double coef = (power / mass) / distance;
        vx += (p->x - this->x) * coef;
        vy += (p->y - this->y) * coef;

    }

    bool isInDoofSkill(std::vector<SkillEffect*>& skillEffects);

    void adjust(std::vector<SkillEffect*>& skillEffects) {
        x = round(x);
        y = round(y);

        if (isInDoofSkill(skillEffects)) {
            vx = round(vx);
            vy = round(vy);
        }
        else {
            vx = round(vx * (1.0 - friction));
            vy = round(vy * (1.0 - friction)); //TODO: remove 1.0 - x
        }
    }

    virtual Collision getCollision() { //why not just return time and already know other stuff
        if (dist(this, &WATERTOWN) + radius >= MAP_RADIUS)
            return Collision(0, this);

        if (vx == 0.0 && vy == 0.0)
            return NULL_COLLISION;

        double a = vx * vx + vy * vy;

        if (a <= 0.0)
            return NULL_COLLISION;

        double b = 2.0 * (x * vx + y * vy);
        double c = x * x + y * y - (MAP_RADIUS - radius) * (MAP_RADIUS - radius);
        double delta = b * b - 4.0 * a * c;

        if (delta <= 0.0)
            return NULL_COLLISION;

        double t = (-b + sqrt(delta)) / (2.0 + a);

        if (t <= 0.0)
            return NULL_COLLISION;

        return Collision(t, this);
    }

    Collision getCollision(Unit* u) {

        if (dist(this, u) <= radius + u->radius)
            return Collision(0.0, this, u);

        if (vx == 0.0 && vy == 0.0 && u->vx == 0.0 && u->vy == 0.0)
            return NULL_COLLISION;

        double x2 = x - u->x;
        double y2 = y - u->y;
        double r2 = radius + u->radius;
        double vx2 = vx - u->vx;
        double vy2 = vy - u->vy;

        double a = vx2 * vx2 + vy2 * vy2;

        if (a <= 0.0)
            return NULL_COLLISION;

        double b = 2.0 * (x2 * vx2 + y2 * vy2);
        double c = x2 * x2 + y2 * y2 - r2 * r2;
        double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
            return NULL_COLLISION;

        double t = (-b - sqrt(delta)) / (2.0 * a);

        if (t <= 0.0)
            return NULL_COLLISION;

        return Collision(t, this, u);
    }

    void bounce(Unit* u) {
        double mcoeff = (mass + u->mass) / (mass * u->mass);
        double nx = x - u->x;
        double ny = y - u->y;
        double nxnysquare = nx * nx + ny * ny;
        double dvx = vx - u->vx;
        double dvy = vy - u->vy;
        double product = (nx * dvx + ny * dvy) / (nxnysquare * mcoeff);
        double fx = nx * product;
        double fy = ny * product;
        double m1c = 1.0 / mass;
        double m2c = 1.0 / u->mass;

        vx -= fx * m1c;
        vy -= fy * m1c;
        u->vx += fx * m2c;
        u->vy += fy * m2c;

        fx = fx * IMPULSE_COEFF;
        fy = fy * IMPULSE_COEFF;

        // Normalize vector at min or max impulse
        double impulse = sqrt(fx * fx + fy * fy);
        double coeff = 1.0;
        if (impulse > EPSILON && impulse < MIN_IMPULSE) {
            coeff = MIN_IMPULSE / impulse;
        }

        fx = fx * coeff;
        fy = fy * coeff;

        vx -= fx * m1c;
        vy -= fy * m1c;
        u->vx += fx * m2c;
        u->vy += fy * m2c;

        double diff = (dist(this, u) - radius - u->radius) / 2.0;
        if (diff <= 0.0) {
            // Unit overlapping. Fix positions.
            moveTo(this, u, diff - EPSILON);
            moveTo(u, this, diff - EPSILON);
        }
    }

    void bounce() {
        double mcoeff = 1.0 / mass;
        double nxnysquare = x * x + y * y;
        double product = (x * vx + y * vy) / (nxnysquare * mcoeff);
        double fx = x * product;
        double fy = y * product;

        vx -= fx * mcoeff;
        vy -= fy * mcoeff;

        fx = fx * IMPULSE_COEFF;
        fy = fy * IMPULSE_COEFF;

        // Normalize vector at min or max impulse
        double impulse = sqrt(fx * fx + fy * fy);
        double coeff = 1.0;
        if (impulse > EPSILON && impulse < MIN_IMPULSE) {
            coeff = MIN_IMPULSE / impulse;
        }

        fx = fx * coeff;
        fy = fy * coeff;
        vx -= fx * mcoeff;
        vy -= fy * mcoeff;

        double diff = dist(this, &WATERTOWN) + radius - MAP_RADIUS;
        if (diff >= 0.0) {
            // Unit still outside of the map, reposition it
            moveTo(this, &WATERTOWN, diff + EPSILON);
        }
    }

};

// ****************************************************************************************

class Tanker : public Unit {
public:
    int water;
    int size;
    Player* player;
    bool killed;

    Tanker(int size, Player* player) : Unit(TYPE_TANKER, 0.0, 0.0) {
        this->player = player;
        this->size = size;
        water = TANKER_EMPTY_WATER;
        mass = TANKER_EMPTY_MASS + TANKER_MASS_BY_WATER * water;
        friction = TANKER_FRICTION;
        radius = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * size;
    }

    Wreck die() {
        if (dist(this, &WATERTOWN) >= MAP_RADIUS)
            return Wreck(0, 0, 0, 0);

        return Wreck(round(x), round(y), water, radius);
    }

    bool isFull() {
        return water >= size;
    }

    void play() {
        if (isFull())
            thrust(&WATERTOWN, -TANKER_THRUST);
        else if (dist(this, &WATERTOWN) > WATERTOWN_RADIUS)
            thrust(&WATERTOWN, TANKER_THRUST);
    }

    Collision getCollision() {
        return NULL_COLLISION;
    }

};
// ****************************************************************************************
class Looter : public Unit {
public:
    int skillCost;
    double skillRange;
    bool skillActive;

    Player* player;

    Point wantedThrustTarget;
    int wantedThrustPower;

    Action attempt;

    Looter(int type, Player* player, double x, double y) : Unit(type, x, y), wantedThrustTarget(0, 0) {
        this->player = player;
        radius = LOOTER_RADIUS;
    }

    SkillEffect* skill(Point* p);

    virtual SkillEffect* skillImpl(Point* p) {
        return nullptr;
    }

    void setWantedThrust(double x, double y, int power) {
        if (power < 0)
            power = 0;
        wantedThrustTarget.x = x;
        wantedThrustTarget.y = y;
        wantedThrustPower = power < MAX_THRUST ? power : MAX_THRUST;
    }
};
// ****************************************************************************************
class Reaper : public Looter {
public:
    Reaper(Player* player, double x, double y) : Looter(LOOTER_REAPER, player, x, y) {
        mass = REAPER_MASS;
        friction = REAPER_FRICTION;
        skillCost = REAPER_SKILL_COST;
        skillRange = REAPER_SKILL_RANGE;
        skillActive = true;
    }

    SkillEffect* skillImpl(Point* p);
};
// ****************************************************************************************
class Destroyer : public Looter {
public:
    Destroyer(Player* player, double x, double y) : Looter(LOOTER_DESTROYER, player, x, y) {
        mass = DESTROYER_MASS;
        friction = DESTROYER_FRICTION;
        skillCost = DESTROYER_SKILL_COST;
        skillRange = DESTROYER_SKILL_RANGE;
        skillActive = true;
    }

    SkillEffect* skillImpl(Point* p);
};
// ****************************************************************************************
class Doof : public Looter {
public:
    Doof(Player* player, double x, double y) : Looter(LOOTER_DOOF, player, x, y) {
        mass = DOOF_MASS;
        friction = DOOF_FRICTION;
        skillCost = DOOF_SKILL_COST;
        skillRange = DOOF_SKILL_RANGE;
        skillActive = true;
    }

    SkillEffect* skillImpl(Point* p);

    int sing() {
        return (int)floor(speed() * DOOF_RAGE_COEF);
    }
};
// ****************************************************************************************
class SkillEffect : public Point {
public:
    int id;
    int type;
    double radius;
    int duration;
    int order;
    bool known;
    Skill skillType;
    Looter* looter;

    SkillEffect(int type, double x, double y, double radius, int duration, int order, Looter* looter) : Point(x, y) {
        id = GLOBAL_ID++;
        this->type = type;
        this->radius = radius;
        this->duration = duration;
        this->looter = looter;
        this->order = order;
    }

    virtual void apply(std::vector<Unit*>& units) {

    }
};

// ****************************************************************************************
class ReaperSkillEffect : public SkillEffect {
public:
    ReaperSkillEffect(int type, double x, double y, double radius, int duration, int order, Reaper* reaper) :
            SkillEffect(type, x , y, radius, duration, order, reaper) {
        skillType = Skill::REAPER;
    }

    void apply(std::vector<Unit*>& units) {
        for (Unit* u : units) {
            if (isInRange(this, u, radius + u->radius))
                u->mass += REAPER_SKILL_MASS_BONUS;
        }
    }
};
// ****************************************************************************************
class DestroyerSkillEffect : public SkillEffect {
public:
    DestroyerSkillEffect(int type, double x, double y, double radius, int duration, int order, Destroyer* destroyer) :
            SkillEffect(type, x, y, radius, duration, order, destroyer) {
        skillType = Skill ::DESTROYER;
    }

    void apply(std::vector<Unit*>& units) {
        for (Unit* u : units) {
            if (isInRange(this, u, radius + u->radius))
                u->thrust(this, -DESTROYER_NITRO_GRENADE_POWER);
        }
    }
};
// ****************************************************************************************
class DoofSkillEffect : public SkillEffect {
public:
    DoofSkillEffect(int type, double x, double y, double radius, int duration, int order, Doof* doof) :
            SkillEffect(type, x, y, radius, duration, order, doof) {
        skillType = Skill ::DOOF;
    }

    void apply(std::vector<Unit*>& units) {
        //No need to do anything
    }
};
// ****************************************************************************************

// FUNCTIONS THAT NEED TO BE DECLARED AFTER

Reaper* Player::getReaper() {
    return dynamic_cast<Reaper*>(looters[LOOTER_REAPER]);
}

Destroyer* Player::getDestroyer() {
    return dynamic_cast<Destroyer*>(looters[LOOTER_DESTROYER]);
}

Doof* Player::getDoof() {
    return dynamic_cast<Doof*>(looters[LOOTER_DOOF]);
}

bool Wreck::harvest(Player* players[], std::vector<SkillEffect*>& skillEffects) {
    for (int i = 0 ; i < 3 ; i++) {
        Player* player = players[i];
        if (isInRange(this, player->getReaper(), radius) && !player->getReaper()->isInDoofSkill(skillEffects)) {
            player->score += 1;
            water -= 1;
        }
    }

    return water > 0;
}

bool Unit::isInDoofSkill(std::vector<SkillEffect*>& skillEffects) {
    for (SkillEffect* skill : skillEffects) {
        if (skill->skillType == DOOF && isInRange(this, skill, skill->radius + radius))
            return true;
    }
    return false;
}

SkillEffect* Looter::skill(Point* p) {
    if (player->rage >= skillCost && dist(this, p) <= skillRange) {
        player->rage -= skillCost;
        return skillImpl(p);
    }
    return nullptr;
}

SkillEffect* Reaper::skillImpl(Point* p) {
    return new ReaperSkillEffect(TYPE_REAPER_SKILL_EFFECT, p->x, p->y, REAPER_SKILL_RADIUS,
                             REAPER_SKILL_DURATION, REAPER_SKILL_ORDER, this);
}

SkillEffect* Destroyer::skillImpl(Point* p) {
    return new DestroyerSkillEffect(TYPE_DESTROYER_SKILL_EFFECT, p->x , p->y,
                                DESTROYER_SKILL_RADIUS, DESTROYER_SKILL_DURATION, DESTROYER_SKILL_ORDER, this);
}

SkillEffect* Doof::skillImpl(Point* p) {
    return new DoofSkillEffect(TYPE_DOOF_SKILL_EFFECT, p->x , p->y, DOOF_SKILL_RADIUS, DOOF_SKILL_DURATION,
                           DOOF_SKILL_ORDER, this);
}

Tanker* Collision::dead() {
    if (a->type == LOOTER_DESTROYER && b->type == TYPE_TANKER && b->mass < REAPER_SKILL_MASS_BONUS) {
        return dynamic_cast<Tanker*>(b);
    }
    if (b->type == LOOTER_DESTROYER && a->type == TYPE_TANKER && a->mass < REAPER_SKILL_MASS_BONUS) {
        return dynamic_cast<Tanker*>(a);
    }
    return nullptr;
}


// ****************************************************************************************
//GLOBAL VARIABLES
Player* players[3];
std::vector<Unit*> units;
std::vector<Looter*> looters;
std::vector<SkillEffect*> skillEffects;
std::vector<Tanker*> tankers;
std::vector<Tanker*> deadTankers;
std::vector<Wreck> wrecks;
// ****************************************************************************************
//GLOBAL METHODS

Collision getNextCollision() {
    Collision result = NULL_COLLISION;

    for (int i = 0; i < units.size(); ++i) {
        Unit* unit = units[i];
        Collision collision = unit->getCollision();

        if (collision.t < result.t)
            result = collision;

        for (int j = i + 1; j < units.size(); ++j) {
            collision = unit->getCollision(units[j]);

            if (collision.t < result.t)
                result = collision;
        }
    }
    return result;
}

void playCollision(Collision& collision) {
    if (collision.b == nullptr) {
        collision.a->bounce();
    }
    else {
        Tanker* dead = collision.dead();

        if (dead != nullptr) {
            deadTankers.push_back(dead);
            tankers.erase(std::remove(tankers.begin(), tankers.end(), dead), tankers.end());
            units.erase(std::remove(units.begin(), units.end(), (Unit*)dead), units.end());

            Wreck wreck = dead->die();

            if (wreck.radius != 0)
                wrecks.push_back(wreck);

        }
    }
}

void updateGame() {
    for (SkillEffect* effect : skillEffects) { //TODO: skillEffects need to be ordered
        effect->apply(units);
    }

    for (Tanker* t : tankers)
        t->play();

    for (Player* player : players) {
        for (Looter* looter : player->looters) {
            if (looter->attempt == Action::MOVE)
                looter->thrust(&looter->wantedThrustTarget, looter->wantedThrustPower);
        }
    }

    double t = 0.0;

    Collision collision = getNextCollision();

    while (collision.t + t <= 1.0) {
        double delta = collision.t;
        for (Unit* unit : units) {
            unit->move(delta);
        }
        t += collision.t;
        playCollision(collision);
        collision = getNextCollision();
    }

    double delta = 1.0 - t;
    for (Unit* unit : units) {
        unit->move(delta);
    }

    //UPDATE TANKER LIST lines 1393-1411 on ref

    for (Unit* unit : units) {
        unit->adjust(skillEffects);
    }

    for (Player* player : players) {
        int increase = player->getDoof()->sing();
        player->rage = MAX_RAGE < player->rage + increase ? MAX_RAGE : player->rage + increase;
    }

    for (Unit* unit : units) {
        while (unit->mass >= REAPER_SKILL_MASS_BONUS)
            unit->mass -= REAPER_SKILL_MASS_BONUS;
    }

    //Remove dead skill effects. Need to design this in a way apply/undo moves work
}
// ****************************************************************************************
int main()
{
    for (int i = 0 ; i < 3; i++) {
        players[i] = new Player(i);
    }
    for (int i = 0 ; i < 3; i++) {
        for (int j = 0 ; j < 3; j++) {
            Looter* looter;
            if (j == LOOTER_REAPER) {
                looter = new Reaper(players[i], 0.0, 0.0);
                units.push_back(looter);
            }
            else if (j == LOOTER_DESTROYER) {
                looter = new Destroyer(players[i], 0.0, 0.0);
                units.push_back(looter);
            }
            else {
                looter = new Doof(players[i], 0.0, 0.0);
                units.push_back(looter);
            }
            looters.push_back(looter);
            players[i]->looters[j] = looter;
        }
    }
    // game loop
    while (true) {
        int myScore;
        std::cin >> myScore; std::cin.ignore();
        int enemyScore1;
        std::cin >> enemyScore1; std::cin.ignore();
        int enemyScore2;
        std::cin >> enemyScore2; std::cin.ignore();
        int myRage;
        std::cin >> myRage; std::cin.ignore();
        int enemyRage1;
        std::cin >> enemyRage1; std::cin.ignore();
        int enemyRage2;
        std::cin >> enemyRage2; std::cin.ignore();
        int unitCount;
        std::cin >> unitCount; std::cin.ignore();
        for (int i = 0; i < unitCount; i++) {
            int unitId;
            int unitType;
            int player;
            float mass;
            int radius;
            int x;
            int y;
            int vx;
            int vy;
            int extra;
            int extra2;
            std::cin >> unitId >> unitType >> player >> mass >> radius >> x >> y >> vx >> vy >> extra >> extra2; std::cin.ignore();
        }

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        std::cout << "WAIT" << std::endl;
        std::cout << "WAIT" << std::endl;
        std::cout << "WAIT" << std::endl;
    }
}
