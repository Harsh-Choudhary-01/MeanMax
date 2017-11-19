#pragma GCC optimize("-O3,inline,omit-frame-pointer,unroll-loops")
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <valarray>
#include <unordered_map>
#include <string>
#include <time.h>
#include <random>
#include <cassert>
#include <assert.h>
#include <chrono>

class Point;
class Unit;
class ActionResult;
class Player;

enum Direction {
    N, NE, E, SE, S, SW, W, NW
};

class Point {
public:
    int x, y, val;

    Point(int x, int y) {
        this->x = x;
        this->y = y;
        this->val = y * 8 + x;
    }

    inline int distance(Point* other) const
    {
        return std::max(abs(x - other->x), abs(y - other->y));
    }

    inline bool operator == (const Point &b) const {
        return b.x == x && b.y == y;
    }

    inline void set(int x1, int y1) {
        x = x1;
        y = y1;
        val = y * 8 + x;
    }
};

struct PointHasher {
    size_t operator()(const Point& obj) const {
        return (31 + obj.x) * 31 + obj.y;
    }
};

struct PointComparator {
    bool operator()(const Point& obj1, const Point& obj2) const {
        return obj1.x == obj2.x && obj1.y == obj2.y;
    }
};

#ifndef WONDEVWOMAN_HEADER_
#define WONDEVWOMEN_HEADER_
void applyMove(ActionResult& result);
void undoMove(ActionResult& result);
int scoreState(bool isPlayer);
std::unordered_set<Point, PointHasher, PointComparator> track(Point& enemyPlace, int prevHeight, bool isZero, Point* end, Point* other);
void getLegalResults(Player& player, std::vector<ActionResult>& results);
int alphabeta(int depth, int alpha, int beta, bool maximizingPlayer, std::chrono::time_point<std::chrono::high_resolution_clock>& startingTime);
int simulate(int depth);
int hash();
#endif

class ActionResult {
public:
    int moveFrom, moveTarget, placeTarget, moveIndex, score, dir1, dir2;
    bool scorePoint, isMove;
    Unit* unit;
    ActionResult(int moveFrom, int moveTarget, int placeTarget, int moveIndex,
                 int dir1, int dir2, bool scorePoint, bool isMove, Unit* unit) : moveFrom(moveFrom), moveTarget(moveTarget), placeTarget(placeTarget),
                                                                                 moveIndex(moveIndex), score(0), dir1(dir1), dir2(dir2), scorePoint(scorePoint), isMove(isMove), unit(unit) {}

    ActionResult() { moveIndex = 5; }

    bool operator < (const ActionResult& res) const
    {
        return (score < res.score);
    }

    bool operator > (const ActionResult& res) const
    {
        return (score > res.score);
    }
};

class Player {
public:
    int score;
    std::vector<Unit> units;
    Player() : score(0), units() {}
};


class Unit {
public:
    int index;
    Player* player;
    Point position;
    bool guessingUnit;

    Unit(Player& player1, int index1) : player(&player1), position(-1, -1) {
        index = index1;
        guessingUnit = false;
    }
};


int arr[64];
int dir[8] = { -8, -7, 1, 9, 8, 7, -1, -9 };
bool occupied[64];
std::string dirValues[8] = { "N", "NE", "E", "SE", "S", "SW", "W", "NW" };
Player playerObj;
Player enemyObj;
int nodeDepthFlood;
int queueArray[64];
int addIndex = 0;
int removeIndex = 0;
int availableAdjacent[4];
int visited[64];
int count = 0;
int MAX_DEPTH = 2;
ActionResult bestRes;
struct timespec startTime, finishTime;
double elapsed;
double totalTime;
const std::string move_string = "MOVE&BUILD";
const std::string push_string = "PUSH&BUILD";
std::unordered_set<Point, PointHasher, PointComparator> possibleEnemyZeroLocations;
std::unordered_set<Point, PointHasher, PointComparator> possibleEnemyOneLocations;
std::vector<int> hashes;
std::unordered_map<long long unsigned, int> hashValues;
int countLeafs = 0;
double timeout = .99;
// Hash Values
long long unsigned table[64][13];
long long unsigned stateHash = 0;
//End Hash Values


inline int fastrand() {
    static unsigned int g_seed = 42;
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

void calculateHash() {
    stateHash = 0;
    for (int i = 0; i < 64; i++) {
        if (!occupied[i]) {
            if (arr[i] != -1)
                stateHash = stateHash ^ table[i][arr[i]];
        }
        else {
            bool enemy = false;
            for (Unit& u : enemyObj.units) {
                if (u.position.val == i) {
                    enemy = true;
                    break;
                }
            }
            if (enemy)
                stateHash = stateHash ^ table[i][9 + arr[i]];
            else
                stateHash = stateHash ^ table[i][5 + arr[i]];
        }
    }
}


int main()
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF);
    for (int i = 0; i < 64; i++)
        for (int j = 0; j < 13; j++)
            table[i][j] = distribution(gen);
    for (int i = 0; i < 64; i++)
        arr[i] = -1;
    for (int i = 0; i < 64; i++)
        occupied[i] = false;
    hashValues.reserve(5000000);
    playerObj.units.emplace_back(playerObj, 0);
    playerObj.units.emplace_back(playerObj, 1);
    int size;
    std::cin >> size;
    std::cin.ignore();
    int unitsPerPlayer;
    std::cin >> unitsPerPlayer;
    std::cin.ignore();
    // game loop
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"
    while (true) {
        countLeafs = 0;
        playerObj.score = 0;
        enemyObj.score = 0;
        bestRes.moveIndex = 5;
        for (Unit& u : enemyObj.units) {
            if (u.guessingUnit)
                occupied[u.position.val] = false;
        }
        enemyObj.units.erase(std::remove_if(enemyObj.units.begin(), enemyObj.units.end(),
                                            [](const Unit& o) { return o.guessingUnit; }), enemyObj.units.end());
        Point* enemyPlace = nullptr;
        int prevHeight = 0;
        for (int i = 0; i < size; i++) {
            std::string row;
            std::cin >> row;
            std::cin.ignore();
            std::cerr << row << "\n";
            for (size_t j = 0; j < row.length(); j++) {
                int ind = i * 8 + j;
                if (row[j] != '.') {
                    int h = row[j] - '0';
                    if (arr[ind] < h && !enemyPlace) {
                        enemyPlace = new Point(j, i);
                        prevHeight = arr[ind];
                    }
                    arr[ind] = h;
                }
            }
        }
        int playerZeroX;
        int playerZeroY;
        int playerOneX;
        int playerOneY;
        int enemyZeroX;
        int enemyZeroY;
        int enemyOneX;
        int enemyOneY;
        std::cin >> playerZeroX >> playerZeroY;
        std::cin.ignore();
        std::cin >> playerOneX >> playerOneY;
        std::cin.ignore();
        std::cin >> enemyZeroX >> enemyZeroY;
        std::cin.ignore();
        std::cin >> enemyOneX >> enemyOneY;
        std::cin.ignore();
        std::cerr << playerZeroX << " " << playerZeroY << "\n";
        std::cerr << playerOneX << " " << playerOneY << "\n";
        std::cerr << enemyZeroX << " " << enemyZeroY << "\n";
        std::cerr << enemyOneX << " " << enemyOneY << "\n";

        if (enemyPlace && (enemyObj.units.size() != 0 || possibleEnemyOneLocations.size() + possibleEnemyZeroLocations.size() != 0) &&
            !(playerObj.units[0].position == *enemyPlace) && !(playerObj.units[1].position == *enemyPlace)) {
            for (Unit& u : enemyObj.units) {
                if (u.index == 0) {
                    possibleEnemyZeroLocations.clear();
                    possibleEnemyZeroLocations.insert(u.position);
                }
                else {
                    possibleEnemyOneLocations.clear();
                    possibleEnemyOneLocations.insert(u.position);
                }
            }
            Point enemyOne(enemyOneX, enemyOneY);
            Point enemyZero(enemyZeroX, enemyZeroY);
            std::unordered_set<Point, PointHasher, PointComparator> list0 = track(*enemyPlace, prevHeight, true,
                                                                                  enemyZeroX == -1 ? nullptr : &enemyZero, enemyOneX == -1 ? nullptr : &enemyOne);
            std::unordered_set<Point, PointHasher, PointComparator> list1 = track(*enemyPlace, prevHeight, false,
                                                                                  enemyOneX == -1 ? nullptr : &enemyOne, enemyZeroX == -1 ? nullptr : &enemyZero);
            if (enemyOneX != -1 && !possibleEnemyOneLocations.count(enemyOne))
                list0.clear();
            if (enemyZeroX != -1 && !possibleEnemyZeroLocations.count(enemyZero))
                list1.clear();
            if (list0.size() == 0) {
                possibleEnemyOneLocations.clear();
                for (const Point& point : list1) {
                    possibleEnemyOneLocations.insert(point);
                }
            }
            else if (list1.size() == 0) {
                possibleEnemyZeroLocations.clear();
                for (const Point& point : list0) {
                    possibleEnemyZeroLocations.insert(point);
                }
            }
            else {
                for (const Point& point : list1) {
                    possibleEnemyOneLocations.insert(point);
                }
                for (const Point& point : list0) {
                    possibleEnemyZeroLocations.insert(point);
                }
            }
        }
        Point playerOne(playerOneX, playerOneY);
        Point playerZero(playerZeroX, playerZeroY);
        for (auto it = std::begin(possibleEnemyZeroLocations); it != std::end(possibleEnemyZeroLocations);) {
            if (it->distance(&playerOne) <= 1 || it->distance(&playerZero) <= 1)
                it = possibleEnemyZeroLocations.erase(it);
            else
                ++it;
        }
        for (auto it = std::begin(possibleEnemyOneLocations); it != std::end(possibleEnemyOneLocations);) {
            if (it->distance(&playerOne) <= 1 || it->distance(&playerZero) <= 1)
                it = possibleEnemyOneLocations.erase(it);
            else
                ++it;
        }
        for (int i = 0; i < 63; i++) {
            occupied[i] = false;
        }
        bool pushed = false;
        int oldPosition = 0;
        int newPosition = 0;
        for (int i = 0; i < unitsPerPlayer; i++) {
            int unitX = i == 0 ? playerZeroX : playerOneX;
            int unitY = i == 0 ? playerZeroY : playerOneY;
            int val = unitY * 8 + unitX;
            if (playerObj.units[i].position.x != -1 && playerObj.units[i].position.val != val) {
                pushed = true;
                oldPosition = playerObj.units[i].position.val;
                newPosition = val;
            }
            playerObj.units[i].position.set(unitX, unitY);
            occupied[val] = true;
        }
        if (!pushed)
            enemyObj.units.clear();
        bool containsOne = false;
        bool containsZero = false;
        for (Unit& u : enemyObj.units) {
            if (u.index == 0)
                containsZero = true;
            else
                containsOne = true;
        }
        for (int i = 0; i < unitsPerPlayer; i++) {
            int otherX = i == 0 ? enemyZeroX : enemyOneX;
            int otherY = i == 0 ? enemyZeroY : enemyOneY;
            if (otherX != -1) {
                occupied[otherY * 8 + otherX] = true;
                if (pushed) {
                    if (i == 0 && containsZero)
                        continue;
                    if (i == 1 && containsOne)
                        continue;
                }
                Unit u(enemyObj, i);
                u.position.set(otherX, otherY);
                enemyObj.units.push_back(u);
            }
        }
        for (Unit& u : enemyObj.units) {
            if (u.index == 0)
                containsZero = true;
            else
                containsOne = true;
        }
        if (!containsZero && possibleEnemyZeroLocations.size() == 1) {
            Unit u(enemyObj, 0);
            u.position.set(possibleEnemyZeroLocations.begin()->x, possibleEnemyZeroLocations.begin()->y);
            enemyObj.units.insert(enemyObj.units.begin(), u);
        }
        if (!containsOne && possibleEnemyOneLocations.size() == 1) {
            Unit u(enemyObj, 1);
            u.position.set(possibleEnemyOneLocations.begin()->x, possibleEnemyOneLocations.begin()->y);
            enemyObj.units.push_back(u);
        }
        if (enemyObj.units.size() != 2 && pushed) {
            int directionTowardOriginOfPush = -100;
            std::unordered_set<Point, PointHasher, PointComparator> list;
            for (int i = 0; i < 8; i++) {
                if (newPosition + dir[i] == oldPosition)
                    directionTowardOriginOfPush = i;
            }
            for (int i = 0; i < 8; i++) {
                int diff = abs(directionTowardOriginOfPush - i);
                if (diff == 1 || diff == 7 || diff == 0) {
                    int newVal = oldPosition + dir[i];
                    list.emplace(newVal % 8, newVal / 8);
                }
            }
            bool oneCanPush = false;
            bool zeroCanPush = false;
            for (const Point& p : possibleEnemyZeroLocations) {
                if (list.count(p))
                    zeroCanPush = true;
            }
            if (containsZero && list.count(enemyObj.units[0].position))
                zeroCanPush = true;
            for (const Point& p : possibleEnemyOneLocations) {
                if (list.count(p))
                    oneCanPush = true;
            }
            if (containsOne && list.count(enemyObj.units[enemyObj.units.size() - 1].position))
                oneCanPush = true;
            if (oneCanPush && !zeroCanPush) {
                for (auto it = std::begin(possibleEnemyOneLocations); it != std::end(possibleEnemyOneLocations);) {
                    if (!list.count(*it))
                        it = possibleEnemyOneLocations.erase(it);
                    else
                        ++it;
                }
                if (!containsOne && possibleEnemyOneLocations.size() == 1) {
                    Unit u(enemyObj, 1);
                    u.position.set(possibleEnemyOneLocations.begin()->x, possibleEnemyOneLocations.begin()->y);
                    enemyObj.units.push_back(u);
                }
            }
            else if (zeroCanPush && !oneCanPush) {
                for (auto it = std::begin(possibleEnemyZeroLocations); it != std::end(possibleEnemyZeroLocations);) {
                    if (!list.count(*it))
                        it = possibleEnemyZeroLocations.erase(it);
                    else
                        ++it;
                }
                if (!containsZero && possibleEnemyZeroLocations.size() == 1) {
                    Unit u(enemyObj, 0);
                    u.position.set(possibleEnemyZeroLocations.begin()->x, possibleEnemyZeroLocations.begin()->y);
                    enemyObj.units.insert(enemyObj.units.begin(), u);
                }
            }
        }
        if (enemyObj.units.size() < 2 && enemyPlace) {
            if (enemyObj.units.size() == 0) {
                std::vector<Point> unitZeroLocations;
                std::vector<Point> unitOneLocations;
                for (const Point& p : possibleEnemyZeroLocations) {
                    unitZeroLocations.push_back(p);
                }
                for (const Point& p : possibleEnemyOneLocations) {
                    unitOneLocations.push_back(p);
                }
                if (unitZeroLocations.size() > 0 && unitOneLocations.size() > 0) {
                    enemyObj.units.emplace_back(enemyObj, 0);
                    enemyObj.units.emplace_back(enemyObj, 1);
                    for (Unit& u : enemyObj.units)
                        u.guessingUnit = true;
                    int minIndex0 = 0;
                    int minIndex1 = 0;
                    int minScore = 2000000;
                    for (size_t i = 0; i < unitZeroLocations.size(); i++) {
                        for (size_t j = 0; j < unitOneLocations.size(); j++) {
                            enemyObj.units[0].position = unitZeroLocations[i];
                            enemyObj.units[1].position = unitOneLocations[j];
                            calculateHash();
                            int score = scoreState(false);
                            if (score < minScore) {
                                minIndex0 = i;
                                minIndex1 = j;
                            }
                        }
                    }
                    enemyObj.units[0].position = unitZeroLocations[minIndex0];
                    enemyObj.units[1].position = unitOneLocations[minIndex1];
                }
                else if (unitZeroLocations.size() > 0 || unitOneLocations.size() > 0) {
                    std::vector<Point>* locations;
                    int index = 0;
                    if (unitZeroLocations.size() > 0)
                        locations = &unitZeroLocations;
                    else {
                        locations = &unitOneLocations;
                        index = 1;
                    }
                    int minIndex = 0;
                    int minScore = 2000000;
                    enemyObj.units.emplace_back(enemyObj, index);
                    enemyObj.units[0].guessingUnit = true;
                    for (size_t i = 0; i < locations->size(); i++) {
                        enemyObj.units[0].position = (*locations)[i];
                        calculateHash();
                        int score = scoreState(false);
                        if (score < minScore) {
                            minScore = score;
                            minIndex = i;
                        }
                    }
                    enemyObj.units[0].position = (*locations)[minIndex];
                }
            }
            else {
                int index = enemyObj.units[0].index == 1 ? 0 : 1;
                if ((index == 0 ? possibleEnemyZeroLocations : possibleEnemyOneLocations).size() > 0) {
                    if (index == 0) {
                        Unit u(enemyObj, index);
                        enemyObj.units.insert(enemyObj.units.begin(), u);
                    }
                    else
                        enemyObj.units.emplace_back(enemyObj, index);
                    enemyObj.units[index].guessingUnit = true;
                    int minScore = 2000000;
                    Point bestPosForEnemy(-1, -1);
                    for (const Point& p : (index == 0 ? possibleEnemyZeroLocations : possibleEnemyOneLocations)) {
                        enemyObj.units[index].position = p;
                        calculateHash();
                        int score = scoreState(false);
                        if (score < minScore) {
                            minScore = score;
                            bestPosForEnemy = p;
                        }
                    }
                    enemyObj.units[index].position = bestPosForEnemy;
                }
            }
        }
        for (Unit& u : enemyObj.units)
            occupied[u.position.val] = true;
        int legalActions;
        std::cin >> legalActions;
        std::cin.ignore();
        std::cerr << legalActions << "\n";
        std::string s;
        for (int i = 0; i < legalActions; i++) {
            std::string atype;
            int index;
            std::string dir1;
            std::string dir2;
            std::cin >> atype >> index >> dir1 >> dir2;
            std::cin.ignore();
            if (i == 0) {
                s.append(atype);
                s.append(" ");
                s.append(std::to_string(index));
                s.append(" ");
                s.append(dir1);
                s.append(" ");
                s.append(dir2);
            }
        }
        stateHash = 0;
        for (int i = 0; i < 64; i++) {
            if (!occupied[i]) {
                if (arr[i] != -1)
                    stateHash = stateHash ^ table[i][arr[i]];
            }
            else {
                bool enemy = false;
                for (Unit& u : enemyObj.units) {
                    if (u.position.val == i) {
                        enemy = true;
                        break;
                    }
                }
                if (enemy)
                    stateHash = stateHash ^ table[i][9 + arr[i]];
                else
                    stateHash = stateHash ^ table[i][5 + arr[i]];
            }
        }
        std::vector<ActionResult> results;
        getLegalResults(playerObj, results);
        MAX_DEPTH = 4;
        nodeDepthFlood = 20;
        if (results.size() < 50)
            nodeDepthFlood = 25;
        if (results.size() < 25)
            nodeDepthFlood = 30;
        if (results.size() <= 30)
            MAX_DEPTH = 5;
        if (results.size() < 10)
            MAX_DEPTH = 6;
        if (timeout > .5)
            MAX_DEPTH = 7;
        results.clear();
        int pos = std::distance(hashes.begin(), find(hashes.begin(), hashes.end(), stateHash));
        if (pos < hashes.size())
            std::cerr << "Collision hash: " << stateHash << " with index: " << pos << "\n";
        std::cerr << "Current state hash " << stateHash << "\n";
        auto startingTime = std::chrono::high_resolution_clock::now();
        if (enemyObj.units.size() > 0) {
            alphabeta(MAX_DEPTH, -2000000, 2000000, true, startingTime);
        }
        else {
            MAX_DEPTH = 2;
            simulate(MAX_DEPTH);
        }
        delete enemyPlace;
        if (enemyObj.units.size()  > 0)
            std::cerr << "Position of enemy 0 " << enemyObj.units[0].position.val << "\n";
        if (enemyObj.units.size() > 1)
            std::cerr << "Position of enemy 1 " << enemyObj.units[1].position.val << "\n";
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - startingTime;
        std::cerr << "Elapsed: " << diff.count() * 1000 << "ms\n";
        totalTime += diff.count() * 1000;
        std::cerr << "Total time: " << totalTime << "\n";
        timeout = .047;
        if (bestRes.moveIndex != 5) {
            std::cerr << "Number leaf nodes: " << countLeafs << "\n";
            applyMove(bestRes);
            std::cout << (bestRes.isMove ? move_string : push_string) << " " << bestRes.moveIndex << " " <<
                      dirValues[bestRes.dir1] << " " << dirValues[bestRes.dir2] << std::endl;
        }
        else {
            std::cout << "ACCEPT-DEFEAT" << std::endl;
        }
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;
    }
#pragma clang diagnostic pop
}

int alphabeta(int depth, int alpha, int beta, bool maximizingPlayer, std::chrono::time_point<std::chrono::high_resolution_clock>& startingTime) {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - startingTime;
    if (diff.count() > timeout)
        return 313131;
    if (depth == 0) {
        countLeafs++;
        return scoreState(maximizingPlayer);
    }
    int size = 6 + 2 * depth;
    bool scoreFlag = false;
    Player& player = maximizingPlayer ? playerObj : enemyObj;
    if (maximizingPlayer) {
        int v = alpha;
        std::vector<ActionResult> results;
        getLegalResults(playerObj, results);
        std::sort(results.begin(), results.end(), [](const ActionResult& lhs, const ActionResult& rhs)
        {
            return lhs.score > rhs.score;
        });
        int s = results.size();
        if (s == 0)
            return scoreState(true) - 2000;
        if (s > size && depth != MAX_DEPTH)
            s = size;
        if (depth == 1) {
            countLeafs += results.size();
            if (results[0].score > v)
                v = results[0].score;
            if (v >= beta)
                return beta;
            return v;
        }
        for (int i = 0; i < s; i++) {
            applyMove(results[i]);
            int vScore = alphabeta(depth - 1, v, beta, false, startingTime);
            if (vScore == 313131) {
                undoMove(results[i]);
                return vScore;
            }
            if (vScore > v) {
                v = vScore;
                if (depth == MAX_DEPTH) {
                    bestRes = results[i];
                }
            }
            undoMove(results[i]);
            if (v >= beta) return beta;
        }
        return v;
    }
    else {
        int v = beta;
        std::vector<ActionResult> results;
        getLegalResults(enemyObj, results);
        std::sort(results.begin(), results.end(), [](const ActionResult& lhs, const ActionResult& rhs)
        {
            return lhs.score < rhs.score;
        });
        int s = results.size();
        if (s == 0)
            return scoreState(false);
        if (s > size)
            s = size;
        if (depth == 1) {
            countLeafs += results.size();
            if (results[0].score < v)
                v = results[0].score;
            if (v <= alpha)
                return alpha;
            return v;
        }
        for (int i = 0; i < s; i++) {
            applyMove(results[i]);
            int vScore = alphabeta(depth - 1, alpha, v, true, startingTime);
            if (vScore == 313131) {
                undoMove(results[i]);
                return vScore;
            }
            if (vScore < v) {
                v = vScore;
            }
            undoMove(results[i]);
            if (v <= alpha) return alpha;
        }
        return v;
    }
}

int simulate(int depth) {
    if (depth == 0)
        return scoreState(true);
    std::vector<ActionResult> results;
    getLegalResults(playerObj, results);
    if (results.size() == 0)
        return scoreState(true) - 2000;
    int maxScore = -20000000;
    for (ActionResult& res : results) {
        applyMove(res);
        int score = simulate(depth - 1);
        if (score > maxScore) {
            if (depth == MAX_DEPTH)
                bestRes = res;
            maxScore = score;
        }
        undoMove(res);
    }
    return maxScore;
}

void computeMove(Unit& unit, int dir1, int dir2, std::vector<ActionResult>& results, bool scoringBool) {

    int positionVal = unit.position.val;
    int targetVal = positionVal + dir[dir1];

    if (targetVal < 0 || targetVal > 63 || arr[targetVal] == -1 || arr[targetVal] >= 4 ||
        arr[targetVal] - arr[positionVal] > 1 || occupied[targetVal])
        return;

    int placeVal = targetVal + dir[dir2];

    if (placeVal < 0 || placeVal > 63 || arr[placeVal] == -1 || arr[placeVal] >= 4 ||
        (occupied[placeVal] && placeVal != positionVal))
        return;

    results.emplace_back(positionVal, targetVal, placeVal, unit.index, dir1, dir2, arr[targetVal] == 3, true, &unit);

    ActionResult& result = results.back();
    applyMove(result);
    result.score = scoreState(scoringBool);
    undoMove(result);
}

void computePush(Unit& unit, int dir1, int dir2, std::vector<ActionResult>& results, bool scoringBool) {
    int diff = abs(dir1 - dir2);

    if (diff != 1 && diff != 7 && diff != 0)
        return;

    int pushVal = unit.position.val + dir[dir1];

    Unit* pushed = nullptr;

    std::vector<Unit>& list = scoringBool ? playerObj.units : enemyObj.units;

    for (Unit& u : list) {
        if (u.position.val == pushVal) {
            pushed = &u;
            break;
        }
    }

    if (!pushed)
        return;

    int resVal = pushVal + dir[dir2];

    if (resVal < 0 || resVal > 63 || occupied[resVal])
        return;

    int toHeight = arr[resVal];

    if (toHeight >= 4 || toHeight == -1 || toHeight - arr[pushVal] > 1)
        return;

    results.emplace_back(pushVal, resVal, pushVal, unit.index, dir1, dir2, false, false, pushed);

    ActionResult& result = results.back();
    applyMove(result);
    result.score = scoreState(scoringBool);
    undoMove(result);
}

void getLegalResults(Player& player, std::vector<ActionResult>& results) {
    for (Unit& unit : player.units) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                computeMove(unit, i, j, results, &player == &enemyObj);
                computePush(unit, i, j, results, &player == &enemyObj);
            }
        }
    }
}

int hash() {
    int h = 64;
    for (int i = 0; i < 63; i++) {
        h = h * 31 + arr[i];
    }
    int a = playerObj.units[0].position.val;
    int b = playerObj.units[1].position.val;
    if (a < b) {
        h = h * 31 + a;
        h = h * 31 + b;
    }
    else {
        h = h * 31 + b;
        h = h * 31 + a;
    }
    if (enemyObj.units.size() == 2) {
        a = enemyObj.units[0].position.val;
        b = enemyObj.units[1].position.val;
        if (a < b) {
            h = h * 31 + a;
            h = h * 31 + b;
        }
        else {
            h = h * 31 + b;
            h = h * 31 + a;
        }
    }
    else if (enemyObj.units.size() == 1) {
        h = h * 31 + enemyObj.units[0].position.val;
    }
    return h;
}

void applyMove(ActionResult& res) {
    occupied[res.unit->position.val] = false;
    stateHash = stateHash ^ table[res.moveFrom][(res.unit->player == &enemyObj ? 9 : 5) + arr[res.moveFrom]];
    stateHash = stateHash ^ table[res.moveTarget][(res.unit->player == &enemyObj ? 9 : 5) + arr[res.moveTarget]];
    stateHash = stateHash ^ table[res.moveTarget][arr[res.moveTarget]];
    res.unit->position.set(res.moveTarget % 8, res.moveTarget / 8);
    if (res.placeTarget != res.moveFrom)
        stateHash = stateHash ^ table[res.placeTarget][arr[res.placeTarget]];
    occupied[res.moveTarget] = true;
    arr[res.placeTarget]++;
    stateHash = stateHash ^ table[res.moveFrom][arr[res.moveFrom]];
    if (res.placeTarget != res.moveFrom)
        stateHash = stateHash ^ table[res.placeTarget][arr[res.placeTarget]];
    if (res.scorePoint)
        res.unit->player->score += 1;
}

void undoMove(ActionResult& res) {
    occupied[res.unit->position.val] = false;
    stateHash = stateHash ^ table[res.moveTarget][(res.unit->player == &enemyObj ? 9 : 5) + arr[res.moveTarget]];
    stateHash = stateHash ^ table[res.moveTarget][arr[res.moveTarget]];
    stateHash = stateHash ^ table[res.moveFrom][arr[res.moveFrom]];
    res.unit->position.set(res.moveFrom % 8, res.moveFrom / 8);
    if (res.placeTarget != res.moveFrom)
        stateHash = stateHash ^ table[res.placeTarget][arr[res.placeTarget]];
    occupied[res.unit->position.val] = true;
    arr[res.placeTarget]--;
    if (res.placeTarget != res.moveFrom)
        stateHash = stateHash ^ table[res.placeTarget][arr[res.placeTarget]];
    stateHash = stateHash ^ table[res.moveFrom][(res.unit->player == &enemyObj ? 9 : 5) + arr[res.moveFrom]];
    if (res.scorePoint)
        res.unit->player->score -= 1;
}

void add(int val) {
    queueArray[addIndex] = val;
    addIndex = (addIndex + 1) > 63 ? 0 : addIndex + 1;
}

int remove() {
    int temp = queueArray[removeIndex];
    removeIndex = (removeIndex + 1) > 63 ? 0 : removeIndex + 1;
    return temp;
}

std::unordered_set<Point, PointHasher, PointComparator> track(Point& enemyPlace, int prevHeight, bool isZero, Point* end, Point* other) {
    arr[enemyPlace.y * 8 + enemyPlace.x] = prevHeight;
    std::unordered_set<Point, PointHasher, PointComparator> list;
    std::unordered_set<Point, PointHasher, PointComparator>& source = isZero ? possibleEnemyZeroLocations : possibleEnemyOneLocations;
    for (const Point& po : source) {
        int p = po.val;
        int height = arr[p];
        for (int change : dir) {
            int loc = p + change;
            int y = loc / 8;
            int x = loc % 8;
            Point location(x, y);
            if (loc > -1 && loc < 64 && arr[loc] != -1 && arr[loc] != 4 && arr[loc] - height <= 1 &&
                location.distance(&enemyPlace) == 1 && !occupied[loc] && (!occupied[enemyPlace.val] || enemyPlace.val == p) &&
                (!end || location == *end) && (!other || !(location == *other))) {
                auto res = list.insert(location);
                std::cerr << "Inserted\n";
            }
        }
    }
    if (source.size() == 0) {
        int p = enemyPlace.val;
        for (int change : dir) {
            int loc = p + change;
            int y = loc / 8;
            int x = loc % 8;
            Point location(x, y);
            if (loc > -1 && loc < 64 && arr[loc] != -1 && arr[loc] != 4 && !occupied[loc] &&
                (!end || location == *end) && (!other || !(location == *other)))
                list.insert(location);
        }
    }
    arr[enemyPlace.y * 8 + enemyPlace.x] = prevHeight + 1;
    return list;
}

int scoreState(bool isPlayer) {
    auto hashVal = hashValues.find(stateHash);
    if (hashVal != hashValues.end())
        return hashVal->second;
    int adjacentNeighbors = 0;
    int scoreDiff = playerObj.score - enemyObj.score;
    for (int i = 0 ; i < 8; i++) {
        for (Unit& u : playerObj.units) {
            int neighborVal = u.position.val + dir[i];
            if (neighborVal < 64 && neighborVal > -1 && arr[neighborVal] != 4 && arr[neighborVal] != -1 && arr[neighborVal] - arr[u.position.val] < 2 && !occupied[neighborVal])
                adjacentNeighbors++;
        }
        for (Unit& u : enemyObj.units) {
            int neighborVal = u.position.val + dir[i];
            if (neighborVal < 64 && neighborVal > -1 && arr[neighborVal] != 4 && arr[neighborVal] != -1 && arr[neighborVal] - arr[u.position.val] < 2 && !occupied[neighborVal])
                adjacentNeighbors--;
        }
    }
    int playerPoints = 0;
    int enemyPoints = 0;
    addIndex = removeIndex = 0;
    for (int i = 0; i < 64; i++)
        visited[i] = -1;
    int p1 = playerObj.units[0].position.val;
    int p2 = playerObj.units[1].position.val;
    bool playerToEnemy = false, enemyToPlayer = false;
    bool isPrevPlayer = false;
    if (isPlayer) {
        add(p1);
        add(p2);
        visited[p1] = 0;
        visited[p2] = 1;
        for (Unit& u : enemyObj.units) {
            add(u.position.val);
            visited[u.position.val] = u.index + 2;
        }
        enemyToPlayer = true;
        isPrevPlayer = true;
    }
    else {
        for (Unit& u : enemyObj.units) {
            add(u.position.val);
            visited[u.position.val] = u.index + 2;
        }
        add(p1);
        add(p2);
        visited[p1] = 0;
        visited[p2] = 1;
        playerToEnemy = true;
    }
    int depth = 0;
    while (addIndex != removeIndex) {
        int p = remove();
        int heightOfPoint = arr[p];
        int index = visited[p];
        bool isPlayerPoint = index < 2;
        if (isPrevPlayer != isPlayerPoint) {
            if (playerToEnemy && !isPlayerPoint && isPrevPlayer)
                depth++;
            else if (enemyToPlayer && isPlayerPoint && !isPrevPlayer)
                depth++;
        }
        if (depth >= nodeDepthFlood)
            break;
        isPrevPlayer = isPlayerPoint;
        for (int i = 0; i < 8; i++) {
            int loc = p + dir[i];
            if (loc > -1 && loc < 64 && arr[loc] != -1 && visited[loc] == -1 && arr[loc] != 4 &&
                arr[loc] - heightOfPoint <= 1) {
                if (isPlayerPoint) {
                    visited[loc] = index;
                    playerPoints++;
                }
                else {
                    visited[loc] = index;
                    enemyPoints++;
                }
                add(loc);
            }
        }
    }
    int s = scoreDiff + adjacentNeighbors + (playerPoints - enemyPoints) * 2;
    hashValues[stateHash] = s;
    return s;
}