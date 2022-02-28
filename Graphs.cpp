#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <stack>

using namespace std;

const int INF = (1 << 30) - 1;
const int MAXXflow = 1001;


class Graph {
private:
    int _n, _m;
    bool _isOriented;
    vector<vector<int>> _adjacentList;

    vector<vector<pair<int, int> >> _adjacentListWithCosts;

    vector<vector<pair<int, int> >> _adjacentListMultigraph;


public:
    Graph(int nodes, int edges, bool oriented) : _n(nodes), _m(edges), _isOriented(oriented) {}

    void addEdges(istream &f); // graf neorientat

    void addEdgesOriented(istream &f); // graf orientat

    void addEdgesWithCosts(ifstream &f); // graf neorientat cu costuri

    void addEdgesOrientedWithCosts(ifstream &f); // graf orientat cu costuri

    vector<int> bfs(int start);

    int nrCompConexe();

    vector<int> sortareTopologica();

    vector<vector<int>> ctc();

    vector<vector<int>> compBiconexe();

    vector<int> Dijkstra(int s);

    vector<pair<int, int>> apm(int node);

    void Bellmanford(int s, ofstream &g);

    void setAdjacentList(const vector<vector<int>> &adjacentList);

    int maxFlow(vector<vector<int>> capacity);

    int darb();

    vector<pair<int, int>> cuplaj();

    int hamilton();

    void addEdgesMultigraph(ifstream &f);

    vector<int> cicluEuler(int start);

private:
    void dfs(int start, vector<int> &visited);

    void dfsCtc(int node, vector<int> &visited, vector<int> &ordine);

    void transpose();

    void biconexe(int node, int parent, stack<int> &s, vector<int> &level, vector<int> &back, vector<vector<int>> &bic);

    bool bfsFlow(int source, int destination, vector<vector<int>> capacity, vector<vector<int>> flow,
                 vector<int> &dad, queue<int> &q);

    void dfsDarb(int start, vector<int> &level, int l, vector<int> &visited);

    bool dfsCuplaj(int node, vector<int> &visited, vector<int> &left, vector<int> &right);
};

void Graph::addEdges(istream &f) {
    int x, y;
    _adjacentList.resize(_n + 1);

    for (int i = 0; i < _m; ++i) {
        f >> x >> y;
        _adjacentList[x].push_back(y);
        _adjacentList[y].push_back(x);
    }
}

void Graph::addEdgesOriented(istream &f) {
    int x, y;
    _adjacentList.resize(_n + 1);

    for (int i = 0; i < _m; ++i) {
        f >> x >> y;
        _adjacentList[x].push_back(y);
    }
}


void Graph::addEdgesOrientedWithCosts(ifstream &f) {
    int x, y, c;
    _adjacentListWithCosts.resize(_n + 1);

    for (int i = 0; i < _m; ++i) {
        f >> x >> y >> c;
        _adjacentListWithCosts[x].push_back(make_pair(y, c));
    }
}

void Graph::addEdgesWithCosts(ifstream &f) {
    int x, y, c;
    _adjacentListWithCosts.resize(_n + 1);

    for (int i = 0; i < _m; ++i) {
        f >> x >> y >> c;
        _adjacentListWithCosts[x].push_back(make_pair(y, c));
        _adjacentListWithCosts[y].push_back(make_pair(x, c));
    }
}

void Graph::addEdgesMultigraph(ifstream &f) {
    int x, y;
    _adjacentListMultigraph.resize(_n + 1);

    for (int i = 1; i <= _m; ++i) {
        f >> x >> y;
        _adjacentListMultigraph[x].push_back(make_pair(y, i));
        _adjacentListMultigraph[y].push_back(make_pair(x, i));
    }
}

void Graph::setAdjacentList(const vector<vector<int>> &adjacentList) {
    _adjacentList = adjacentList;
}

// parametri: nodul de start
// output: vector ce contine distanta din nodul sursa start la celelalte noduri
// Se insereaza nodul start intr-o coada vida, cu costul 0. La fiecare pas, se ia nodul din inceputul cozii,
// se elimina si apoi se adauga vecinii nevizitati la finalul cozii. Costul unui nod nou adaugat va fi costul nodului
// care l-a adaugat + 1
vector<int> Graph::bfs(int start) {
    queue<int> q;
    vector<int> dist(_n + 1);

    for (int i = 1; i <= _n; i++) {
        dist[i] = -1;
    }

    q.push(start);
    dist[start] = 0;
    while (!q.empty()) {
        int node = q.front();

        for (int i = 0; i < _adjacentList[node].size(); i++) {
            if (dist[_adjacentList[node][i]] == -1) {
                q.push(_adjacentList[node][i]);
                dist[_adjacentList[node][i]] = dist[node] + 1;
            }
        }
        q.pop();
    }

    return dist;
}

// parametri: nodul, vectorul de vizitari
void Graph::dfs(int node, vector<int> &visited) {
    visited[node] = 1;
    for (int i = 0; i < _adjacentList[node].size(); i++) {
        if (!visited[_adjacentList[node][i]])
            dfs(_adjacentList[node][i], visited);
    }
}

// output: numarul componentelor conexe
int Graph::nrCompConexe() {
    int i, nr = 0;
    vector<int> visited(_n + 1, 0);

    for (i = 1; i <= _n; i++) {
        if (!visited[i]) {
            nr++;
            dfs(i, visited);
        }
    }
    return nr;
}


// parametri: fisierul de intrare
// output: true daca se poate forma un graf cu secventa de grade data, false altfel
bool havel_hakimi(istream &f) {
    vector<int> grade;
    int x, s = 0;

    //citesc gradele din fisier si calculez si suma
    while (f >> x) {
        grade.push_back(x);
        s += x;
    }

    //verific sa nu fie un nod cu gradul mai mare de n-1
    for (int i = 0; i < grade.size(); i++)
        if (grade[i] > grade.size() - 1) {
            return false;
        }

    //verific sa nu fie suma impara
    if (s % 2) {
        return false;
    }

    //sortez gradele descrescator
    sort(grade.begin(), grade.end(), greater<int>());

    //cat timp gradele nodurilor >0, rulam algoritmul
    while (grade[0] != 0) {
        //scad cu 1 gradele urmatoarelor grad[0] noduri, apoi grad[0] devine 0
        for (int i = 1; i <= grade[0]; i++) {
            grade[i]--;
            if (grade[i] < 0) {
                return false;
            }
        }
        grade[0] = 0;
        sort(grade.begin(), grade.end(), greater<int>()); // la fiecare pas reordonam gradele descrescator
    }

    return true;
}


// output: vector ce contine sortarea topologica a nodurilor grafului dat
vector<int> Graph::sortareTopologica() {
    vector<int> grade_interne(_n + 2, 0);
    queue<int> coada;
    vector<int> answer;

    //calculez gradele interne ale nodurilor
    for (int i = 1; i <= _n; i++)
        for (int j = 0; j < _adjacentList[i].size(); j++)
            grade_interne[_adjacentList[i][j]] += 1;

    //adaug in coada toate nodurile cu gradul intern 0
    for (int i = 1; i <= _n; i++)
        if (grade_interne[i] == 0)
            coada.push(i);

    int varf;
    while (!coada.empty()) {
        //extragem un varf din coada
        varf = coada.front();
        answer.push_back(varf);
        coada.pop();

        //scad gradele interne ale vecinilor, fara sa elimin nodul din reprezentare si adaug in coada vecinii al caror
        // grad intern devine 0
        for (int i = 0; i < _adjacentList[varf].size(); i++) {
            grade_interne[_adjacentList[varf][i]]--;
            if (grade_interne[_adjacentList[varf][i]] == 0)
                coada.push(_adjacentList[varf][i]);
        }
    }

    return answer;
}

//graful transpus
void Graph::transpose() {
    vector<int> transpEdges[_n + 1];
    int i, j;
    for (i = 1; i <= _n; i++)
        for (j = 0; j < _adjacentList[i].size(); j++)
            transpEdges[_adjacentList[i][j]].push_back(i);
    for (i = 1; i <= _n; i++) {
        _adjacentList[i] = transpEdges[i];
    }
}

// parametri: nodul, vectorul de vizitari, vectorul de ordine
void Graph::dfsCtc(int node, vector<int> &visited, vector<int> &ordine) {
    visited[node] = 1;
    for (int i = 0; i < _adjacentList[node].size(); i++) {
        if (!visited[_adjacentList[node][i]])
            dfsCtc(_adjacentList[node][i], visited, ordine);
    }
    //in ordine pastrez nodurile cand sunt finalizate de dfs
    ordine.push_back(node);
}

// output: matrice care contine pe fiecare linie o componentă tare conexă prin enumerarea nodurilor componente
vector<vector<int>> Graph::ctc() {
    vector<vector<int>> ans;
    vector<int> visited(_n + 1, 0);

    vector<int> ordine;
    for (int i = 1; i <= _n; i++) {
        if (!visited[i]) {
            dfsCtc(i, visited, ordine);
        }
    }
    //inversez ordinea
    reverse(ordine.begin(), ordine.end());

    //fac graful transpus
    transpose();

    //resetez vectorul de vizitari
    for (int i = 1; i <= _n; i++)
        visited[i] = 0;

    //fac dfs pentru elemente din ordine (cele inversate)
    for (int i = 0; i < _n; i++) {
        if (!visited[ordine[i]]) {
            vector<int> v;
            dfsCtc(ordine[i], visited, v);
            ans.push_back(v);
        }
    }

    return ans;
}

// parametri: nodul, parintele nodului, stiva, vector pentru nivelul fiecarui nod, vector pentru minimul de muchii
// prin care ne putem intoarce fara a trece prin parinte, matricea care contine pe fiecare linie cate o componenta biconexa
void
Graph::biconexe(int node, int parent, stack<int> &s, vector<int> &level, vector<int> &back, vector<vector<int>> &bic) {
    //calculez nivelurile nodurilor
    level[node] = level[parent] + 1;
    back[node] = level[node];

    s.push(node);
    //parcurg vecinii
    for (int i = 0; i < _adjacentList[node].size(); i++) {
        if (level[_adjacentList[node][i]]) {
            if (_adjacentList[node][i] != parent)
                //actualizez minimul
                back[node] = min(level[_adjacentList[node][i]], back[node]);
        } else {
            biconexe(_adjacentList[node][i], node, s, level, back, bic);  //dfs pentru biconexe
            back[node] = min(back[_adjacentList[node][i]], back[node]);
            if (level[node] <= back[_adjacentList[node][i]]) { //verific daca face parte din componenta biconexa
                vector<int> v;

                while (s.top() != _adjacentList[node][i]) {
                    v.push_back(s.top()); //pun in v elementele componentei biconexe
                    s.pop();
                }

                v.push_back(_adjacentList[node][i]);
                s.pop();
                v.push_back(node);

                bic.push_back(v);
            }
        }
    }
}

// output: matrice care contine pe fiecare linie cate o componenta biconexa
vector<vector<int>> Graph::compBiconexe() {
    vector<vector<int>> bic;

    //initializez cu 0 vectorii pentru nivel si pentru intoarcere
    vector<int> level(_n + 1, 0);
    vector<int> back(_n + 1, 0);
    stack<int> s;

    for (int i = 1; i <= _n; i++) {
        //pentru fiecare nod, in functia biconexe incep sa calculez nivelul daca acesta nu este inca setat
        if (!level[i]) {
            biconexe(i, 0, s, level, back, bic);
            //golesc stiva
            while (!s.empty()) {
                s.pop();
            }
        }
    }

    return bic;
}

// parametri: nodul de start
// output:  vector cu n-1 numere naturale separate printr-un spatiu. Al i-lea numar va reprezenta lungimea unui
// drum minim de la nodul start la nodul i+1
vector<int> Graph::Dijkstra(int start) {
    vector<int> dist(_n + 1, INF);
    vector<int> visited(_n + 1, 0);

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;

    dist[start] = 0;
    q.push(make_pair(0, start));

    while (!q.empty()) {
        int currentNode = q.top().second; //extrag varful cu eticheta minima din q
        q.pop();

        if (visited[currentNode] == 0) {
            visited[currentNode] = 1;

            for (int i = 0; i < _adjacentListWithCosts[currentNode].size(); ++i) {
                int neighbour = _adjacentListWithCosts[currentNode][i].first;
                int cost = _adjacentListWithCosts[currentNode][i].second;

                if (dist[neighbour] > dist[currentNode] + cost) { // actualizez drumul cand gasesc un cost mai mic
                    dist[neighbour] = dist[currentNode] + cost;
                    q.push(make_pair(dist[neighbour], neighbour));
                }
            }
        }

    }
    return dist;
}

// parametri: nodul de start
// output: vector de pair care contine muchiile ce apartin arborelui solutie
vector<pair<int, int>> Graph::apm(int start) {
    vector<pair<int, int>> answer;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
    vector<int> costs(_n + 1, INF);
    vector<int> dad(_n + 1, 0);
    vector<int> visited(_n + 1, 0);
    int total = 0;
    dad[start] = 0;
    costs[start] = 0;
    q.push(make_pair(0, start));

    while (!q.empty()) {
        int currentNode = q.top().second;
        q.pop();

        if (visited[currentNode] == 0) {
            visited[currentNode] = 1;

            for (int i = 0; i < _adjacentListWithCosts[currentNode].size(); ++i) {
                int neighbour = _adjacentListWithCosts[currentNode][i].first;
                int cost = _adjacentListWithCosts[currentNode][i].second;
                if (visited[neighbour] == 0 && cost < costs[neighbour]) { // actualizez cand gasesc un cost mai mic la un vecin nevizitat
                    costs[neighbour] = cost;
                    dad[neighbour] = currentNode;
                    q.push(make_pair(costs[neighbour], neighbour));
                }
            }
            total += costs[currentNode];
            if (dad[currentNode] != 0)
                answer.push_back(make_pair(currentNode, dad[currentNode]));
        }
    }

    return answer;
}

// parametri: nodul x, vectorul de tati
// output: radacina nodului
int Find(int x, vector<int> &dad) {
    int root, aux;
    root = x;

    while (root != dad[root])
        root = dad[root];

    while (x != dad[x]) {
        aux = dad[x];
        dad[x] = root;
        x = aux;
    }

    return root;
}

// parametri: nodul x, nodul y, vectorul de tati, vectorul de ranguri
// output: realizeaza reuniunea dintre multimile in care se afla elementul x, respectiv elementul y
void Union(int x, int y, vector<int> &dad, vector<int> &rang) {
    int root_x = Find(x, dad);
    int root_y = Find(y, dad);

    if (rang[root_x] > rang[root_y]) {
        dad[root_y] = root_x;
        rang[root_x] += rang[root_y];
    } else {
        dad[root_x] = root_y;
        rang[root_y] += rang[root_x];
    }

}

// parametri: nodul de start, fisierul in care se scrie rezultatul
// output: in fisier este scris costul minim al unui lanţ de la nodul start la celelalte noduri sau "ciclu negativ" daca apare
void Graph::Bellmanford(int start, ofstream &g) {
    vector<int> dist(_n + 1, INF);
    queue<int> q;  // optimizare folosind coada
    vector<int> inQueue(_n + 1, 0);
    vector<int> visited(_n + 1, 0);
    dist[start] = 0;
    q.push(start);
    inQueue[start] = 1;

    while (!q.empty()) {
        int currentNode = q.front();
        q.pop();
        visited[currentNode] += 1; // numar de cate ori a fost vizitat un nod pentru a detecta ciclurile
        inQueue[currentNode] = 0;

        if (visited[currentNode] >= _n) {
            g << "Ciclu negativ!";
            return;
        }


        for (int i = 0; i < _adjacentListWithCosts[currentNode].size(); ++i) {
            int neighbour = _adjacentListWithCosts[currentNode][i].first;
            int cost = _adjacentListWithCosts[currentNode][i].second;

            if (dist[neighbour] > dist[currentNode] + cost) { // actualizez cand gasesc un cost mai mic
                dist[neighbour] = dist[currentNode] + cost;
                if (inQueue[neighbour] == 0) q.push(neighbour);
            }
        }

    }

    for (int i = 2; i <= _n; ++i)
        g << dist[i] << " ";

}

// parametri: nodul sursa, nodul destinatie, matricea pentru capacitati, matricea pentru flow, vector de tati, coada
// output: true daca s-a ajuns la destinatie, false altfel
bool Graph::bfsFlow(int source, int destination, vector<vector<int>> capacity, vector<vector<int>> flow,
                    vector<int> &dad, queue<int> &q) {
    for (int i = 1; i <= _n; i++) {
        dad[i] = 0;
    }

    q.push(source);
    dad[source] = source;
    int arrivedAtDestination = 0;

    while (!q.empty()) {
        int currentNode = q.front();
        q.pop();

        if (currentNode == destination) {
            arrivedAtDestination = 1;
            continue;
        }

        for (auto &node: _adjacentList[currentNode]) {
            if ((capacity[currentNode][node] - flow[currentNode][node] > 0) && dad[node] == 0) {
                dad[node] = currentNode;
                q.push(node);
            }
        }
    }
    return arrivedAtDestination;
}

// parametri: matrice pentru capacitati
// output: capacitatea maxima
int Graph::maxFlow(vector<vector<int>> capacity) {
    int source = 1, destination = _n;
    vector<vector<int>> flow(MAXXflow, vector<int>(MAXXflow, 0));
    vector<int> dad(_n + 1, 0);
    queue<int> q;
    int answer = 0;

    while (bfsFlow(source, destination, capacity, flow, dad, q)) {
        for (auto &node: _adjacentList[destination]) {
            int currentNode = destination;
            int minim = INF;

            if (dad[node] == 0) {
                continue;
            }

            dad[destination] = node;
            while (currentNode != source) {
                minim = min(minim, capacity[dad[currentNode]][currentNode] - flow[dad[currentNode]][currentNode]);
                currentNode = dad[currentNode];
            }

            currentNode = destination;
            while (currentNode != source) {
                flow[dad[currentNode]][currentNode] += minim;
                flow[currentNode][dad[currentNode]] -= minim;
                currentNode = dad[currentNode];
            }
            answer += minim;
        }
    }

    return answer;
}

// parametri: nodul, vector pentru nivel, vector de vizitari
void Graph::dfsDarb(int node, vector<int> &level, int l, vector<int> &visited) {
    visited[node] = 1;
    level[node] = l;
    for (int i = 0; i < _adjacentList[node].size(); ++i) {
        if (!visited[_adjacentList[node][i]]) {
            dfsDarb(_adjacentList[node][i], level, l + 1, visited);
        }
    }
}

// output: diametrul arborelui
int Graph::darb() {
    int diameter = 0;
    vector<int> level(_n + 1, 0);
    vector<int> visited(_n + 1, 0);
    int maxLevel = 0;
    int finale = 0;

    dfsDarb(1, level, 1, visited);

    for (int i = 0; i < _n; ++i)
        if (level[i] > maxLevel) {
            maxLevel = level[i];
            finale = i;
        }

    for (int i = 0; i <= _n; ++i) {
        visited[i] = 0;
    }
    level.clear();
    dfsDarb(finale, level, 1, visited);

    for (int i = 0; i <= _n; ++i)
        if (level[i] > diameter)
            diameter = level[i];

    return diameter;
}

// parametri: numarul de noduri, matricea in care se construiesc drumurile minime
// output: se formeaza matricea drumurilor minime
void royfloyd(int n, int matrix[105][105]) {

    for (int k = 1; k <= n; ++k) {
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                if (matrix[i][j] > matrix[i][k] + matrix[k][j])
                    matrix[i][j] = matrix[i][k] + matrix[k][j];
    }
}

// parametri: nodul, vectorul de vizitari, perechile nodurilor din stanga, perechile nodurilor din dreapta
// output: true, daca s-a format pereche sau lant alternant, si false altfel
// DFS-ul acesta va cauta nodurile din prima multime care inca nu au pereche si o va forma
// sau nodurile care au pereche din prima multime, dar continua cu lant alternant care accepta o augmentare de o unitate de flux
// Functia returneaza true, daca s-a format pereche sau lant alternant, si false altfel.
bool Graph::dfsCuplaj(int node, vector<int> &visited, vector<int> &left, vector<int> &right) {

    if (visited[node])
        return false;

    visited[node] = 1;

    for (int i = 0; i < _adjacentList[node].size(); ++i)
        if (!right[_adjacentList[node][i]] || dfsCuplaj(right[_adjacentList[node][i]], visited, left, right)) {
            right[_adjacentList[node][i]] = node;
            left[node] = _adjacentList[node][i];
            return true;
        }

    return false;
}

// output: vector de pair care reprezinta cuplajul
vector<pair<int, int>> Graph::cuplaj() {
    vector<pair<int, int>> answer;

    vector<int> left(10001, 0), right(10001, 0);
    vector<int> visited(_n + 1, 0);
    bool ok;

    do {
        ok = false;
        for (int i = 0; i <= _n; ++i)
            visited[i] = 0;

        for (int i = 1; i <= _n; ++i)
            if (left[i] == 0 && dfsCuplaj(i, visited, left, right))
                ok = true;
    } while (ok == true);

    for (int i = 1; i <= _n; ++i)
        if (left[i] != 0)
            answer.push_back(make_pair(i, left[i]));

    return answer;
}

// output: costul minim al nodurilor care ajung in nodul 0
// Se creeaza matricea costurilor a tuturor drumurilor, care initial e initializata cu INF.
// Se foloseste reprezentarea binara pentru a retine nodurile care fac parte din lant, acestea fiind marcate cu 1
// Pentru fiecare lant, se verifica nodurile care fac parte din acesta si se actualizeaza costul minim
// La final este cautat costul minim al nodurilor care ajung in nodul 0, deoarece cu acesta am inceput
int Graph::hamilton() {
    int answer = INF;
    int nr_nodes = 1 << _n; //numarul de noduri
    int cost[nr_nodes][_n];

    for (int i = 0; i < nr_nodes; ++i)
        for (int j = 0; j < _n; ++j)
            cost[i][j] = INF;

    cost[1][0] = 0;  // fixez nodul de inceput 0, care are costul 0

    for (int i = 0; i < nr_nodes; ++i)
        for (int j = 0; j < _n; ++j)
            if ((i & (1 << j))) {
                for (int k = 0; k < _adjacentListWithCosts[j].size(); ++k) {

                    if (i & (1 << _adjacentListWithCosts[j][k].first)) {
                        cost[i][j] = min(cost[i][j], cost[i ^ (1 << j)][_adjacentListWithCosts[j][k].first] +
                                                     _adjacentListWithCosts[j][k].second);
                    }
                }
            }

    for (int i = 0; i < _adjacentListWithCosts[0].size(); ++i)
        answer = min(answer, cost[nr_nodes - 1][_adjacentListWithCosts[0][i].first] +
                             _adjacentListWithCosts[0][i].second);

    return answer;
}


// parametri: nodul de start
// output: vector ce contine ciclul eulerian
// Se verifica intai daca toate nodurile au gradele pare.
// La inceput se pune intr-o stiva nodul de start.
// Cat timp stiva nu este vida, se ia nodul curent, apoi ultimul vecin al acestuia din lista,
// pentru a fi usor de sters muchia.
// Se sterge muchia cu pop_back(), apoi se viziteaza muchia si se adauga in stiva nodul la care a dus muchia.
// Daca nodul nu mai are vecini, atunci este scos din stiva si adaugat la ciclul eulerian.
vector<int> Graph::cicluEuler(int start) {
    vector<int> answer;
    vector<int> visitedEdges(_m + 1, 0);

    for (int i = 0; i < _n; ++i)
        if (_adjacentListMultigraph[i].size() % 2 != 0) {
            answer.push_back(-1);
            return answer;
        }

    stack<int> euler;
    euler.push(start);

    while (!euler.empty()) {
        int currentNode = euler.top();

        if (_adjacentListMultigraph[currentNode].size() != 0) {
            int neighbour = _adjacentListMultigraph[currentNode].back().first;
            int index = _adjacentListMultigraph[currentNode].back().second;
            _adjacentListMultigraph[currentNode].pop_back();

            if (visitedEdges[index] == 0) {
                visitedEdges[index] = 1;
                euler.push(neighbour);
            }
        } else {
            euler.pop();
            answer.push_back(currentNode);
        }
    }

    return answer;
}


void solve_bfs() {
    ifstream f("bfs.in");
    ofstream g("bfs.out");

    int n, m, s;
    f >> n >> m >> s;

    Graph gr(n, m, true);
    gr.addEdgesOriented(f);

    vector<int> answer = gr.bfs(s);
    for (int i = 1; i <= n; ++i) {
        g << answer[i] << " ";
    }

    f.close();
    g.close();
}

void solve_nrCompConexe() {
    ifstream f("dfs.in");
    ofstream g("dfs.out");
    int n, m;
    f >> n >> m;
    Graph gr(n, m, false);
    gr.addEdges(f);

    g << gr.nrCompConexe();

    f.close();
    g.close();
}

void solve_sortareTopologica() {
    ifstream f("sortaret.in");
    ofstream g("sortaret.out");
    int n, m;
    f >> n >> m;
    Graph gr(n, m, true);
    gr.addEdgesOriented(f);

    vector<int> answer = gr.sortareTopologica();

    for (int i = 0; i < answer.size(); ++i)
        g << answer[i] << ' ';

    f.close();
    g.close();
}

void solve_ctc() {
    ifstream f("ctc.in");
    ofstream g("ctc.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, true);
    gr.addEdgesOriented(f);

    vector<vector<int>> answer = gr.ctc();

    g << answer.size() << '\n';

    for (int i = 0; i < answer.size(); i++) {
        for (int j = 0; j < answer[i].size(); j++)
            g << answer[i][j] << " ";
        g << '\n';
    }
}

void solve_biconex() {
    ifstream f("biconex.in");
    ofstream g("biconex.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, false);
    gr.addEdges(f);

    vector<vector<int>> answer = gr.compBiconexe();

    g << answer.size() << '\n';

    for (int i = 0; i < answer.size(); i++) {
        for (int j = 0; j < answer[i].size(); j++)
            g << answer[i][j] << " ";
        g << '\n';
    }

    f.close();
    g.close();
}

void solve_Dijkstra() {
    ifstream f("dijkstra.in");
    ofstream g("dijkstra.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, true);
    gr.addEdgesOrientedWithCosts(f);

    vector<int> answer = gr.Dijkstra(1);

    for (int i = 2; i <= n; ++i)
        if (answer[i] != INF)
            g << answer[i] << ' ';
        else g << 0 << ' ';
}

void solve_disjoint() {
    ifstream f("disjoint.in");
    ofstream g("disjoint.out");

    int n, m, cod, x, y;
    f >> n >> m;

    vector<int> dad(n + 1);
    vector<int> rang(n + 1, 0);

    for (int i = 1; i <= n; ++i) {
        dad[i] = i;
    }

    for (int i = 0; i < m; ++i) {
        f >> cod >> x >> y;
        if (cod == 1) Union(x, y, dad, rang);
        else {
            if (Find(x, dad) == Find(y, dad))
                g << "DA\n";
            else g << "NU\n";
        }
    }

    f.close();
    g.close();
}

void solve_HavelHakimi() {
    ifstream f("hh.in");
    ofstream g("hh.out");

    bool answer = havel_hakimi(f);

    if (answer) g << "Da.\n";
    else g << "Nu.\n";
    f.close();
}

void solve_apm() {
    ifstream f("apm.in");
    ofstream g("apm.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, false);
    gr.addEdgesWithCosts(f);
    vector<pair<int, int>> answer = gr.apm(1);

    g << answer.size() << '\n';
    for (int i = 0; i < answer.size(); ++i) {
        g << answer[i].first << " " << answer[i].second << '\n';
    }

    f.close();
    g.close();
}

void solve_bellmanford() {
    ifstream f("bellmanford.in");
    ofstream g("bellmanford.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, false);
    gr.addEdgesOrientedWithCosts(f);

    gr.Bellmanford(1, g);

    f.close();
    g.close();
}

void solve_maxFlow() {
    ifstream f("maxflow.in");
    ofstream g("maxflow.out");

    int n, m, answer, x, y, z;
    f >> n >> m;
    Graph grf(n, m, false);
    vector<vector<int>> gr(100001);
    vector<vector<int>> capacity(MAXXflow, vector<int>(MAXXflow, 0));

    for (int i = 0; i < m; ++i) {
        f >> x >> y >> z;
        gr[x].push_back(y);
        gr[y].push_back(x);
        capacity[x][y] = z;
    }
    grf.setAdjacentList(gr);

    answer = grf.maxFlow(capacity);
    g << answer;

    f.close();
    g.close();
}

void solve_darb() {
    ifstream f("darb.in");
    ofstream g("darb.out");

    int n, x, y;
    f >> n;
    Graph gr(n, n - 1, false);
    gr.addEdges(f);

    g << gr.darb();

    f.close();
    g.close();
}

void solve_royfloyd() {
    ifstream f("royfloyd.in");
    ofstream g("royfloyd.out");

    int n;
    int matrix[105][105];
    f >> n;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            f >> matrix[i][j];
            if (matrix[i][j] == 0 && i != j)
                matrix[i][j] = 1001;
        }
    }

    royfloyd(n, matrix);

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j)
            if (matrix[i][j] == 1001)
                g << 0 << " ";
            else
                g << matrix[i][j] << "  ";
        g << '\n';
    }

    f.close();
    g.close();

}

void solve_cicluEuler() {
    ifstream f("ciclueuler.in");
    ofstream g("ciclueuler.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, false);
    gr.addEdgesMultigraph(f);
    vector<int> answer = gr.cicluEuler(1);
    if (answer[0] == -1) {
        g << -1;
    } else {
        for (int i = 0; i < answer.size() - 1; ++i) {
            g << answer[i] << " ";
        }
    }
    f.close();
    g.close();
}

void solve_hamilton() {
    ifstream f("hamilton.in");
    ofstream g("hamilton.out");

    int n, m;
    f >> n >> m;
    Graph gr(n, m, true);
    gr.addEdgesOrientedWithCosts(f);
    int answer = gr.hamilton();
    if (answer == INF)
        g << "Nu exista solutie";
    else
        g << answer;

    f.close();
    g.close();
}

void solve_cuplaj() {
    ifstream f("cuplaj.in");
    ofstream g("cuplaj.out");

    int n, m, e, x, y;
    f >> n >> m >> e;

    Graph gr(n, e, true);

    gr.addEdgesOriented(f);

    vector<pair<int, int>> answer = gr.cuplaj();

    g << answer.size() << '\n';
    for (int i = 0; i < answer.size(); ++i)
        g << answer[i].first << " " << answer[i].second << '\n';

    f.close();
    g.close();
}

int main() {
    solve_HavelHakimi();
    solve_royfloyd();
    solve_disjoint();
    solve_cicluEuler();
    solve_Dijkstra();
    solve_apm();
    solve_bellmanford();
    solve_hamilton();
    solve_bfs();
    solve_nrCompConexe();
    solve_sortareTopologica();
    solve_biconex();
    solve_ctc();
    solve_cuplaj();
    solve_darb();
    solve_maxFlow();
    return 0;
}
