#pragma GCC optimize("-O3,inline,omit-frame-pointer,unroll-loops")
#include <iostream>
#include <chrono>
#include <string>
#include <algorithm>
#include <math.h>
#include <vector>
#include <unordered_set>
#include <set>
#include <fstream>

using namespace std::chrono;

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

constexpr int DEPTH = 4; //TODO: refactor to put DEPTH in usage
constexpr int POOL = 50;
constexpr double MUTATION = 2;

// Global first free id for all elements on the map 
int GLOBAL_ID = 0;
int temp_ID = 0;

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

high_resolution_clock::time_point start;
#define NOW high_resolution_clock::now()
#define TIME duration_cast<duration<double>>(NOW - start).count()

// ***********************************************************

static unsigned int g_seed;
inline void fast_srand(int seed) {
    //Seed the generator
    g_seed = seed;
}
inline int fastrand() {
    //fastrand routine returns one integer, similar output value range as C lib.
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}
inline int fastRandInt(int maxSize) {
    return fastrand() % maxSize;
}
inline int fastRandInt(int a, int b) {
    return(a + fastRandInt(b - a));
}
inline double fastRandDouble() {
    return static_cast<double>(fastrand()) / 0x7FFF;
}
inline double fastRandDouble(double a, double b) {
    return a + (static_cast<double>(fastrand()) / 0x7FFF)*(b-a);
}

// ****************************************************************************************
//SAVED VARIABLES
std::vector<SkillEffect*> savedSkillEffects;
std::vector<SkillEffect*> skillsToDelete;
std::vector<Tanker*> savedTankers;
std::vector<Wreck*> savedWrecks;
std::vector<Wreck*> wrecksToDelete;
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

    double dx = u2->x - u1->x;
    double dy = u2->y - u1->y;
    double coef = distance / d;

    u1->x += dx * coef;
    u1->y += dy * coef;
}

inline bool isInRange(Point* p1, Point* p2, double range) {
    return p1 != p2 && dist(p1, p2) <= range;
}
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
    int savedDuration;

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

    void save() {
        savedDuration = duration;
    }

    void reset() {
        duration = savedDuration;
    }
};
// ****************************************************************************************
struct SkillEffectComparator {
    bool operator()(SkillEffect* obj1, SkillEffect* obj2) const {
        int order = obj1->order - obj2->order;
        if (order != 0)
            return order < 0;
        return obj1->id - obj2->id < 0;
    }
};
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

};

// ****************************************************************************************

// Center of the map
Point WATERTOWN(0, 0);

Collision NULL_COLLISION(1.0 + EPSILON);
double NULL_COLLISION_TIME = 1.0 + EPSILON;

// ****************************************************************************************
class Player {
public:
    int score;
    int index;
    int rage;
    Looter* looters[3];
    bool dead;
    int saveValues[2];

    Player(int index) {
        this->index = index;
    }

    void kill() {
        dead = true;
    }

    void save();

    void reset();

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
    int savedWater;

    Wreck(double x, double y, int water, double radius) : Point(x, y) {
        id = GLOBAL_ID++;
        this->radius = radius;
        this->water = water;
    }

    bool harvest(Player* players[], std::set<SkillEffect*, SkillEffectComparator>& skillEffects);

    void save() {
        savedWater = water;
    }

    void reset() {
        water = savedWater;
    }

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
    double saveValues[5];

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

    virtual void save() {
        saveValues[0] = x;
        saveValues[1] = y;
        saveValues[2] = vx;
        saveValues[3] = vy;
        saveValues[4] = mass;
    }

    virtual void reset() {
        x = saveValues[0];
        y = saveValues[1];
        vx = saveValues[2];
        vy = saveValues[3];
        mass = saveValues[4];
    }

    void thrust(Point* p, int power) {
        double distance = dist(this, p);

        if (distance <= EPSILON)
            return;

        double coef = (power / mass) / distance;
        vx += (p->x - this->x) * coef;
        vy += (p->y - this->y) * coef;

    }

    bool isInDoofSkill(std::set<SkillEffect*, SkillEffectComparator>& skillEffects);

    void adjust(std::set<SkillEffect*, SkillEffectComparator>& skillEffects) {
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

    virtual double getCollision() { //why not just return time and already know other stuff
        if (dist(this, &WATERTOWN) + radius >= MAP_RADIUS)
            return 0;

        if (vx == 0.0 && vy == 0.0)
            return NULL_COLLISION_TIME;

        double a = vx * vx + vy * vy;

        if (a <= 0.0)
            return NULL_COLLISION_TIME;

        double b = 2.0 * (x * vx + y * vy);
        double c = x * x + y * y - (MAP_RADIUS - radius) * (MAP_RADIUS - radius);
        double delta = b * b - 4.0 * a * c;

        if (delta <= 0.0)
            return NULL_COLLISION_TIME;

        double t = (-b + sqrt(delta)) / (2.0 + a);

        if (t <= 0.0)
            return NULL_COLLISION_TIME;

        return t;
    }

    double getCollision(Unit* u) {

        if (dist(this, u) <= radius + u->radius)
            return 0;

        if (vx == 0.0 && vy == 0.0 && u->vx == 0.0 && u->vy == 0.0)
            return NULL_COLLISION_TIME;

        double x2 = x - u->x;
        double y2 = y - u->y;
        double r2 = radius + u->radius;
        double vx2 = vx - u->vx;
        double vy2 = vy - u->vy;

        double a = vx2 * vx2 + vy2 * vy2;

        if (a <= 0.0)
            return NULL_COLLISION_TIME;

        double b = 2.0 * (x2 * vx2 + y2 * vy2);
        double c = x2 * x2 + y2 * y2 - r2 * r2;
        double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
            return NULL_COLLISION_TIME;

        double t = (-b - sqrt(delta)) / (2.0 * a);

        bool test = t < EPSILON;

        if (t <= 0.0)
            return NULL_COLLISION_TIME;

        return t;
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
    int savedWater;

    Tanker(int size, Player* player) : Unit(TYPE_TANKER, 0.0, 0.0) {
        this->player = player;
        this->size = size;
        water = TANKER_EMPTY_WATER;
        mass = TANKER_EMPTY_MASS + TANKER_MASS_BY_WATER * water;
        friction = TANKER_FRICTION;
        radius = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * size;
    }

    Wreck* die() {
        if (dist(this, &WATERTOWN) >= MAP_RADIUS)
            return nullptr;

        Wreck* temp = new Wreck(round(x), round(y), water, radius);
        wrecksToDelete.push_back(temp);
        return temp;
    }

    void save() {
        Unit::save();
        savedWater = water;
    }

    void reset() {
        Unit::reset();
        water = savedWater;
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

    double getCollision() {
        return NULL_COLLISION_TIME;
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
        attempt = Action::MOVE;
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
class ReaperSkillEffect : public SkillEffect {
public:
    ReaperSkillEffect(int type, double x, double y, double radius, int duration, int order, Reaper* reaper) :
            SkillEffect(type, x , y, radius, duration, order, reaper) {
        skillType = Skill::REAPER;
    }

    void apply(std::vector<Unit*>& units) {
        duration -= 1;
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
        duration -= 1;
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
        duration -= 1;
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

bool Wreck::harvest(Player* players[], std::set<SkillEffect*, SkillEffectComparator>& skillEffects) {
    for (int i = 0 ; i < 3 ; i++) {
        Player* player = players[i];
        if (isInRange(this, player->getReaper(), radius) && !player->getReaper()->isInDoofSkill(skillEffects)) {
            player->score += 1;
            water -= 1;
        }
    }

    return water > 0;
}

bool Unit::isInDoofSkill(std::set<SkillEffect*, SkillEffectComparator>& skillEffects) {
    for (SkillEffect* skill : skillEffects) {
        if (skill->skillType == DOOF && isInRange(this, skill, skill->radius + radius))
            return true;
    }
    return false;
}

SkillEffect* Looter::skill(Point* p) {
    if (player->rage >= skillCost && dist(this, p) <= skillRange) {
        player->rage -= skillCost;
        SkillEffect* effect = skillImpl(p);
        skillsToDelete.push_back(effect);
        return effect;
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

void Player::save() {
    saveValues[0] = score;
    saveValues[1] = rage;
    for (Looter* looter : looters)
        looter->save();
}

void Player::reset() {
    score = saveValues[0];
    rage =saveValues[1];
    for (Looter* looter : looters)
        looter->reset();
}

struct SkillEffectHasher {
    size_t operator()(const SkillEffect& obj) const {
        return (31 + obj.id);
    }
};

struct TankerHasher {
    size_t operator()(const Tanker& obj) const {
        return (31 + obj.id);
    }
};

struct TankerComparator {
    bool operator()(const Tanker& obj1, const Tanker& obj2) const {
        return obj1.id == obj2.id;
    }
};

struct UnitHasher {
    size_t operator()(const Unit* obj) const {
        return (31 + obj->id);
    }
};

struct UnitComparator {
    bool operator()(const Unit* obj1, const Unit* obj2) const {
        return obj1->id == obj2->id;
    }
};

struct WreckHasher {
    size_t operator()(const Wreck& obj) const {
        return (31 + obj.id);
    }
};

struct WreckComparator {
    bool operator()(const Wreck& obj1, const Wreck& obj2) const {
        return obj1.id == obj2.id;
    }
};


// ****************************************************************************************
//GLOBAL VARIABLES
Player* players[3];
std::vector<Unit*> units;
std::vector<Looter*> looters;
std::set<SkillEffect*, SkillEffectComparator> skillEffects;
std::vector<Tanker*> tankers;
std::vector<Wreck*> wrecks;
int turn = 0;
Unit* aObject = nullptr;
Unit* bObject = nullptr;
// ****************************************************************************************
//GLOBAL METHODS

Tanker* dead() {
    if (aObject->type == LOOTER_DESTROYER && bObject->type == TYPE_TANKER && bObject->mass < REAPER_SKILL_MASS_BONUS) {
        return dynamic_cast<Tanker*>(bObject);
    }
    if (bObject->type == LOOTER_DESTROYER && aObject->type == TYPE_TANKER && aObject->mass < REAPER_SKILL_MASS_BONUS) {
        return dynamic_cast<Tanker*>(aObject);
    }
    return nullptr;
}

double getNextCollision() {
    double result = NULL_COLLISION_TIME;

    for (int i = 0; i < units.size(); ++i) {
        Unit* unit = units[i];
        double collision = unit->getCollision();

        if (collision < result) {
            result = collision;
            aObject = unit;
            bObject = nullptr;
        }

        for (int j = i + 1; j < units.size(); ++j) {
            collision = unit->getCollision(units[j]);

            if (collision < result) {
                result = collision;
                aObject = unit;
                bObject = units[j];
            }
        }
    }
    return result;
}

void playCollision() {
    if (bObject == nullptr) {
        aObject->bounce();
    }
    else {
        Tanker* deadTanker = dead();

        if (deadTanker != nullptr) {
            tankers.erase(std::remove(tankers.begin(), tankers.end(), deadTanker), tankers.end());
            units.erase(std::remove(units.begin(), units.end(), (Unit*)deadTanker), units.end());

            Wreck* wreck = deadTanker->die();

            if (wreck != nullptr)
                wrecks.push_back(wreck);

        }
        else {
            aObject->bounce(bObject);
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

    double collision = getNextCollision();

    while (collision + t <= 1.0) {
        double delta = collision;
        for (Unit* unit : units) {
            unit->move(delta);
        }
        t += collision;
        playCollision();
        collision = getNextCollision();
    }

    double delta = 1.0 - t;
    for (Unit* unit : units) {
        unit->move(delta);
    }

    //UPDATE TANKER LIST lines 1393-1411 on ref
    std::unordered_set<Unit*, UnitHasher, UnitComparator> tankersToRemove;

    for (Tanker* tanker : tankers) {
        double distance = dist(tanker, &WATERTOWN);
        bool full = tanker->isFull();

        if (distance <= WATERTOWN_RADIUS && !full) {
            tanker->water += 1;
            tanker->mass += TANKER_MASS_BY_WATER;
        }
        else if (distance >= TANKER_SPAWN_RADIUS + tanker->radius && full)
            tankersToRemove.insert(tanker);
    }

    units.erase(std::remove_if(units.begin(), units.end(),
                               [&tankersToRemove](Unit* o) { return tankersToRemove.count(o); }), units.end());

    tankers.erase(std::remove_if(tankers.begin(), tankers.end(),
                                 [&tankersToRemove](Tanker* o) { return tankersToRemove.count(o); }), tankers.end());

    wrecks.erase(std::remove_if(wrecks.begin(), wrecks.end(),
                                [](Wreck* o) { bool alive = o->harvest(players, skillEffects);
                                    return !alive;}), wrecks.end());

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

    for (auto it = skillEffects.begin(); it != skillEffects.end();) {
        if ((*it)->duration <= 0)
            it = skillEffects.erase(it);
        else
            ++it;
    }

    //Remove dead skill effects. Need to design this in a way apply/undo moves work
}


void save() {
    skillsToDelete.clear();
    wrecksToDelete.clear();
    savedSkillEffects.clear();
    savedTankers.clear();
    savedWrecks.clear();

    temp_ID = GLOBAL_ID;

    for (SkillEffect* effect : skillEffects) {
        savedSkillEffects.push_back(effect);
        effect->save();
    }

    for (Tanker* tanker : tankers) {
        savedTankers.push_back(tanker);
        tanker->save();
    }

    for (Wreck* wreck : wrecks) {
        savedWrecks.push_back(wreck);
        wreck->save();
    }

    for (Player* player : players)
        player->save();

}

void reset() {

    units.resize(9);

    skillEffects.clear();
    tankers.clear();
    wrecks.clear();

    for (SkillEffect* effect : skillsToDelete)
        delete effect;

    for (Wreck* wreck : wrecksToDelete)
        delete wreck;

    for (SkillEffect* effect : savedSkillEffects) {
        skillEffects.insert(effect);
        effect->reset();
    }

    for (Tanker* tanker : savedTankers) {
        tankers.push_back(tanker);
        tanker->reset();
    }

    for (Wreck* wreck : savedWrecks) {
        wrecks.push_back(wreck);
        wreck->reset();
    }

    for (Player* player : players)
        player->reset();

    skillsToDelete.clear();
    wrecksToDelete.clear();

    GLOBAL_ID = temp_ID;
}


void heuristic(Player* player) { //Sets player moves based on heuristic
    double dist1 = 50000;
    double dist2 = 0;

    Wreck* closest = nullptr;
    for (Wreck* wreck : wrecks) {
        dist2 = dist(player->looters[0], wreck);
        if (dist2 < dist1) {
            closest = wreck;
            dist1 = dist2;
        }
    }

    if (closest == nullptr)
        player->looters[0]->setWantedThrust(player->looters[1]->x, player->looters[1]->y, 300);
    else
        player->looters[0]->setWantedThrust(closest->x, closest->y, dist1 < 800 ? 0 : 300);

    dist1 = 50000;

    Tanker* closestTanker = nullptr;
    for (Tanker* tanker: tankers) {
        dist2 = dist(player->looters[1], tanker);
        if (dist2 < dist1) {
            closestTanker = tanker;
            dist1 = dist2;
        }
    }

    if (closestTanker == nullptr)
        player->looters[1]->setWantedThrust(0, 0, 300);
    else
        player->looters[1]->setWantedThrust(closestTanker->x, closestTanker->y, 300);

    player->looters[2]->attempt = Action::WAIT;
}

int scoreState() {
    int score = players[0]->score - players[1]->score - players[2]->score;
    score *= 1000;
    double dist1 = dist(players[0]->looters[2], players[1]->looters[0]);
    double dist2 = dist(players[0]->looters[2], players[2]->looters[0]);
    score -= (int)(dist1 < dist2 ? dist1 : dist2);
    dist1 = 50000;
    for (Wreck* wreck : wrecks) {
        dist2 = dist(players[0]->looters[0], wreck);
        dist1 = (dist1 < dist2 ? dist1 : dist2);
    }
    if (dist1 == 50000)
        dist1 = dist(players[0]->looters[0], players[0]->looters[1]);
    score -= (int)dist1;
    dist1 = 50000;
    for (Tanker* tanker : tankers) {
        dist2 = dist(players[0]->looters[1], tanker);
        dist1 = (dist1 < dist2 ? dist1 : dist2);
    }
    if (dist1 == 50000)
        dist1 = dist(players[0]->looters[1], 0, 0);
    score -= (int)dist1;
    return score;
}

void print() {
    for (int i = 0 ; i < 3 ; i++)
        std::cout << players[i]->score << std::endl;
    for (int i = 0 ; i < 3; i++)
        std::cout << players[i]->rage << std::endl;
    for (Unit* unit : units)
        std::cout << unit->id << " " << unit->type << " " << unit->mass
                  << " " << unit->radius << " " << unit->x << " " << unit->y
                  << " " << unit->vx << " " << unit->vy << std::endl;

}

class Move {
public:
    int moveType;
    int x;
    int y;
    int thrust;

    Move() {

    }

    void randomize() {
        x = fastRandInt(-120, 120);
        y = fastRandInt(-120, 120);
        thrust = fastRandInt(0, 500);
    }

    void mutate(double amplitude) {
        //X
        double minAmp = x - 20 * amplitude;
        double maxAmp = x + 20 * amplitude;
        if (minAmp < -120)
            minAmp = -120;
        if (maxAmp > 120)
            maxAmp = 120;
        x = fastRandInt(minAmp, maxAmp);

        //Y
        minAmp = y - 20 * amplitude;
        maxAmp = y + 20 * amplitude;
        if (minAmp < -120)
            minAmp = -120;
        if (maxAmp > 120)
            maxAmp = 120;
        y = fastRandInt(minAmp, maxAmp);

        minAmp = thrust - 25 * amplitude;
        maxAmp = thrust + 25 * amplitude;
        if (minAmp < 0)
            minAmp = 0;
        if (maxAmp > 500)
            maxAmp = 500;
        thrust = fastRandInt(minAmp, maxAmp);
    }

};


class Solution {
public:
    Move movesReaper[4];
    Move movesDestroyer[4];
    Move movesDoof[4];
    int score;

    Solution() {

    }

    void randomize() {
        for (Move& move : movesReaper)
            move.randomize();
        for (Move& move : movesDestroyer)
            move.randomize();
        for (Move& move : movesDoof)
            move.randomize();
    }

    Solution* mutate(double amplitude) {
        Solution* solution = copy();
        for (Move& move : solution->movesReaper)
            move.mutate(amplitude);
        for (Move& move : solution->movesDestroyer)
            move.mutate(amplitude);
        for (Move& move : solution->movesDoof)
            move.mutate(amplitude);
        return solution;
    }

    void mutateInPlace(double amplitude) {
        for (Move& move : movesReaper)
            move.mutate(amplitude);
        for (Move& move : movesDestroyer)
            move.mutate(amplitude);
        for (Move& move : movesDoof)
            move.mutate(amplitude);
    }

    Solution* merge(Solution* other) {
        Solution* child = new Solution();
        for (int i = 0 ; i < 4; i++) {
            if (fastRandInt(2))
                child->movesReaper[i] = movesReaper[i];
            else
                child->movesReaper[i] = other->movesReaper[i];
            if (fastRandInt(2))
                child->movesDoof[i] = movesDoof[i];
            else
                child->movesDoof[i] = other->movesDoof[i];
            if (fastRandInt(2))
                child->movesDestroyer[i] = other->movesDestroyer[i];
            else
                child->movesDestroyer[i] = movesDestroyer[i];
        }
        return child;
    }

    void mergeInPlace(Solution* other) {
        for (int i = 0 ; i < 4; i++) {
            if (fastRandInt(2)) {
                movesReaper[i] = other->movesReaper[i];
                movesDestroyer[i] = other->movesDestroyer[i];
                movesDoof[i] = other->movesDoof[i];
            }
        }
    }

    Solution* copy() {
        Solution* copy = new Solution();
        for (int i = 0 ; i < 4; i++) {
            copy->movesReaper[i].x = movesReaper[i].x;
            copy->movesReaper[i].y = movesReaper[i].y;
            copy->movesReaper[i].thrust = movesReaper[i].thrust;
            //copy->movesReaper[i].moveType = movesReaper[i].moveType;

            copy->movesDestroyer[i].x = movesDestroyer[i].x;
            copy->movesDestroyer[i].y = movesDestroyer[i].y;
            copy->movesDestroyer[i].thrust = movesDestroyer[i].thrust;
            //copy->movesDestroyer[i].moveType = movesDestroyer[i].moveType;

            copy->movesDoof[i].x = movesDoof[i].x;
            copy->movesDoof[i].y = movesDoof[i].y;
            copy->movesDoof[i].thrust = movesDoof[i].thrust;
            //copy->movesDoof[i].moveType = movesDoof[i].moveType;
        }
        copy->score = score;
        return copy;
    }

    void copy(Solution* other) {
        for (int i = 0 ; i < 4; i++) {
            movesReaper[i].x = other->movesReaper[i].x;
            movesReaper[i].y = other->movesReaper[i].y;
            movesReaper[i].thrust = other->movesReaper[i].thrust;
            //movesReaper[i].moveType = other->movesReaper[i].moveType;

            movesDestroyer[i].x = other->movesDestroyer[i].x;
            movesDestroyer[i].y = other->movesDestroyer[i].y;
            movesDestroyer[i].thrust = other->movesDestroyer[i].thrust;
            //movesDestroyer[i].moveType = other->movesDestroyer[i].moveType;

            movesDoof[i].x = other->movesDoof[i].x;
            movesDoof[i].y = other->movesDoof[i].y;
            movesDoof[i].thrust = other->movesDoof[i].thrust;
            //movesDoof[i].moveType = other->movesDoof[i].moveType;
        }
        score = other->score;
    }


    void simulate() {
        for (int i = 0; i < 4; i++) {
            players[0]->looters[0]->setWantedThrust(movesReaper[i].x * 50,
                                                    movesReaper[i].y * 50, movesReaper[i].thrust);
            players[0]->looters[1]->setWantedThrust(movesDestroyer[i].x * 50,
                                                    movesDestroyer[i].y * 50, movesDestroyer[i].thrust);
            players[0]->looters[2]->setWantedThrust(movesDoof[i].x * 50,
                                                    movesDoof[i].y * 50, movesDoof[i].thrust);
            heuristic(players[1]);
            heuristic(players[2]);
            updateGame();
        }
        score = scoreState();
        reset();
    }

};
// ****************************************************************************************
int main()
{
    //std::ifstream in("~/CLionProjects/MeanMax/input.txt");
    //std::cin.rdbuf(in.rdbuf());
    fast_srand(42);
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

    Solution* best = new Solution();
    // game loop
    while (true) {
        // ****************************************************************************************
        //Resetting state
        GLOBAL_ID = 0;
        for (SkillEffect* skillEffect : skillEffects)
            delete skillEffect;
        for (Tanker* tanker : tankers)
            delete tanker;
        for (Wreck* wreck : wrecks)
            delete wreck;
        skillEffects.clear();
        tankers.clear();
        wrecks.clear();
        units.resize(9);
        // ****************************************************************************************
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

        ///*
        std::cerr << myScore << std::endl;
        std::cerr << enemyScore1 << std::endl;
        std::cerr << enemyScore2 << std::endl;
        std::cerr << myRage << std::endl;
        std::cerr << enemyRage1 << std::endl;
        std::cerr << enemyRage2 << std::endl;
        std::cerr << unitCount << std::endl;
        //*/

        players[0]->score = myScore;
        players[0]->rage = myRage;
        players[1]->score = enemyScore1;
        players[1]->rage = enemyRage1;
        players[2]->score = enemyScore2;
        players[2]->rage = enemyRage2;
        int tempID = 0;
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
            ///*
            std::cerr << unitId << " " << unitType << " " << player << " " << mass << " " << radius << " " << x << " " <<
                      y << " " << vx << " " << vy << " " << extra << " " << extra2 << std::endl;
            //*/
            tempID = unitId < tempID ? tempID : unitId;
            if (unitType < 3) {
                players[player]->looters[unitType]->mass = mass;
                players[player]->looters[unitType]->x = x;
                players[player]->looters[unitType]->y = y;
                players[player]->looters[unitType]->vx = vx;
                players[player]->looters[unitType]->vy = vy;
                players[player]->looters[unitType]->radius = radius;
                players[player]->looters[unitType]->id = unitId;
            }
            else if (unitType == 3) {
                Tanker* tanker = new Tanker(extra2, players[player]);
                tanker->mass = mass;
                tanker->water = extra;
                tanker->size = extra2;
                tanker->radius = radius;
                tanker->x = x;
                tanker->y = y;
                tanker->vx = vx;
                tanker->vy = vy;
                tanker->id = unitId;
                tankers.push_back(tanker);
                units.push_back(tanker);
            }
            else if (unitType == 4) {
                Wreck* wreck = new Wreck(x, y, extra, radius);
                wreck->id = unitId;
                wrecks.push_back(wreck);
            }
            else if (unitType == 5) { //JUST CREATE SKILLS THROUGH YOUR PLAYER. DOESN'T MATTER

                SkillEffect* skillEffect = new ReaperSkillEffect(TYPE_REAPER_SKILL_EFFECT, x, y,
                                                                 REAPER_SKILL_RADIUS, REAPER_SKILL_DURATION,
                                                                 REAPER_SKILL_ORDER, players[0]->getReaper());
                skillEffect->duration = extra;
                skillEffect->id = unitId;
                skillEffects.insert(skillEffect);
            }
            else {
                SkillEffect* skillEffect = new DoofSkillEffect(TYPE_DOOF_SKILL_EFFECT, x , y, DOOF_SKILL_RADIUS,
                                                               DOOF_SKILL_DURATION, DOOF_SKILL_ORDER, players[0]->getDoof());
                skillEffect->duration = extra;
                skillEffect->id = unitId;
                skillEffects.insert(skillEffect);
            }
        }
        GLOBAL_ID = tempID + 1;

        save(); //SAVING STATE

        start = NOW;

        //TODO: use previous GA solution and modify the last turn randomly

        Solution* base;

        if (turn) {
            base = new Solution();

            for (int j = 1; j < 4; ++j) {
                base->movesReaper[j - 1] = best->movesReaper[j];
                base->movesDestroyer[j - 1] = best->movesDestroyer[j];
                base->movesDoof[j - 1] = best->movesDoof[j];
            }

            delete best;
        }

        Solution** pool = new Solution*[POOL];
        Solution** newPool = new Solution*[POOL];
        Solution** temp;
        int counter = POOL;

        best = new Solution();
        Solution* sol = new Solution();
        sol->randomize();

        sol->simulate();
        pool[0] = sol;

        best->copy(sol);

        Solution* tempBest = sol;

        int startI = 1;

        //TODO: Fill 1/5 pool with mutations of last turn if not turn 0

        if (turn) {
            for (int i = startI; i < POOL / 5; ++i) {
                Solution* solution = new Solution();
                solution->copy(base);

                solution->movesReaper[3].randomize();
                solution->movesDestroyer[3].randomize();
                solution->movesDoof[3].randomize();

                solution->simulate();

                if (solution->score > tempBest->score)
                    tempBest = solution;

                pool[i] = solution;
            }

            delete base;

            startI = POOL / 5;
        }

        for (int i = startI; i < POOL; ++i) { //Fill rest with totally random
            Solution* solution = new Solution();
            solution->randomize();

            solution->simulate();

            if (solution->score > tempBest->score)
                tempBest = solution;

            pool[i] = solution;
        }

        if (tempBest->score > best->score)
            best->copy(tempBest);

        tempBest = best;

        double limit = turn ? .043 : .9;

#define LIMIT TIME < limit

        bool continueLoop = true;

        int poolFE;
        while (LIMIT) {

            Solution* solution = new Solution();
            solution->copy(tempBest);
            solution->mutateInPlace(.05);
            solution->simulate();

            if (solution->score > tempBest->score)
                tempBest = solution;

            newPool[0] = solution;

            counter += 1;

            poolFE = 1;

            while (poolFE < POOL && LIMIT) {
                int aIndex = fastRandInt(POOL);
                int bIndex;

                do {
                    bIndex = fastRandInt(POOL);
                } while (bIndex == aIndex);

                int firstIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                do {
                    aIndex = fastRandInt(POOL);
                } while (aIndex == firstIndex);

                do {
                    bIndex = fastRandInt(POOL);
                } while (bIndex == aIndex || bIndex == firstIndex);

                int secondIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                Solution* child = pool[firstIndex]->merge(pool[secondIndex]);

                //if (!fastRandInt(MUTATION))
                //child->mutateInPlace(.2);

                child->simulate();

                if (child->score > tempBest->score)
                    tempBest = child;

                newPool[poolFE++] = child;

                counter += 1;
            }

            for (int i = 0 ; i < POOL; ++i)
                delete pool[i];

            temp = pool;
            pool = newPool;
            newPool = temp;

            if (tempBest->score > best->score) {
                best->copy(tempBest);
            }

            tempBest = best;

        }

        reset();

        std::cout << best->movesReaper[0].x * 50 << " " << best->movesReaper[0].y * 50 << " " <<
                  best->movesReaper[0].thrust << std::endl;
        std::cout << best->movesDestroyer[0].x * 50 << " " << best->movesDestroyer[0].y * 50 <<
                  " " << best->movesDestroyer[0].thrust << std::endl;
        std::cout << best->movesDoof[0].x * 50 << " " << best->movesDoof[0].y * 50 << " "
                  << best->movesDoof[0].thrust << std::endl;

        std::cerr << "Counter: " << counter << std::endl;

        for (int i = 0; i < poolFE; ++i)
            delete pool[i];

        delete [] pool;
        delete [] newPool;


        turn += 1;


        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

//        std::cout << "0 0 300" << std::endl;
//        if (players[0]->rage >= 60)
//            std::cout << "SKILL 0 0" << std::endl;
//        else
//            std::cout << "0 0 300" << std::endl;
//        std::cout << "0 0 300" << std::endl;
    }
}
