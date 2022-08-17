#include "Grafo.h"
#include "Vertice.h"
#include "Aresta.h"
#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <list>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <float.h>
#include <iomanip>
#include <string.h>
#include <vector>
#include <algorithm>
#include <map>
#include <chrono>

using namespace std;

/**************************************************************************************************
 * Defining the Graph's methods
 **************************************************************************************************/

// Constructor
Grafo::Grafo(int order, bool directed, bool weighted_edge, bool weighted_node) {
    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_node = weighted_node;
    this->first_node = this->last_node = nullptr;
    this->number_edges = 0;
}

Grafo::Grafo(int order, int cluster, string clusterType, vector<tuple<int, int>> clustersLimits) {
    // this->order = order;
    this->cluster = cluster;
    this->clusterType = clusterType;
    this->clustersLimits = clustersLimits;
    this->first_node = this->last_node = nullptr;
    this->weighted_edge = 1;
    this->weighted_node = 1;
}

Grafo::Grafo(int order, int cluster, double clustersCapacity) {
    // this->order = order;
    this->cluster = cluster;
    this->clustersCapacity = clustersCapacity;
    this->first_node = this->last_node = nullptr;
    this->weighted_edge = 1;
    this->weighted_node = 1;
}


Grafo::Grafo(int inferiorLimit, int upperLimit) {
    this->first_node = this->last_node = nullptr;
    this->weighted_edge = 1;
    this->weighted_node = 1;
    this->inferiorLimit = inferiorLimit;
    this->upperLimit = upperLimit;
    this->currentLimit = 0;
}
// Destructor
Grafo::~Grafo() {
    Vertice *next_node = this->first_node;
    while (next_node != nullptr) {
        next_node->removeAllEdges();
        Vertice *aux_node = next_node->getNextNode();
        delete next_node;
        next_node = aux_node;
    }
}
long inline size_arq(fstream &arq) {
    arq.ignore(std::numeric_limits<std::streamsize>::max());
    std::streamsize size= arq.gcount();
    arq.clear();
    arq.seekg(0, std::ios_base::beg);

    return size;
}

// Getters
int Grafo::getOrder() {
    return this->order;
}

int Grafo::getNumberEdges() {
    return this->number_edges;
}

// Function that verifies if the graph is directed
bool Grafo::getDirected() {
    return this->directed;
}

// Function that verifies if the graph is weighted at the edges
bool Grafo::getWeightedEdge() {
    return this->weighted_edge;
}

// Function that verifies if the graph is weighted at the nodes
bool Grafo::getWeightedNode() {
    return this->weighted_node;
}

Vertice *Grafo::getFirstNode() {
    return this->first_node;
}

Vertice *Grafo::getLastNode() {
    return this->last_node;
}


void Grafo::insertNode(int id) {
    Vertice *next;
    Vertice *aux = nullptr;

    // Verifica se já existe algum nó
    if (this->getFirstNode() == nullptr) {
        this->first_node = new Vertice(id);
        this->last_node = this->getFirstNode();
    } else {
        if (!this->searchNode(id)) {
            next = this->getFirstNode();
            // Procura o último nó inserido
            while (next != nullptr) {
                aux = next;
                next = next->getNextNode();
            }
            // Inseri o nó na última posição
            aux->setNextNode(new Vertice(id));
            this->last_node = this->getNode(id);
        }
    }
}

void Grafo::insertNodeAndWeight(int id, int weight) {
    Vertice *next;
    Vertice *aux = nullptr;

    // Verifica se já existe algum nó
    if (this->getFirstNode() == nullptr) {
        this->first_node = new Vertice(id);
        // cout << "144 - id:" << this->first_node->getId() << endl;
        this->first_node->setWeight(weight);
        // this->last_node = this->getFirstNode();
        this->order++;
    } else {
        if (!this->searchNode(id)) {
            next = this->getFirstNode();
            // Procura o último nó inserido
            while (next != nullptr) {
                aux = next;
                next = next->getNextNode();
            }
            // Inseri o nó na última posição
            aux->setNextNode(new Vertice(id));
            this->last_node = this->getNode(id);
            this->last_node->setWeight(weight);
            // cout << "TEST" << endl;
            this->order++;
        }
    }
}

void Grafo::insertEdge(int id, int target_id, float weight) {
    // Procura se o nó id existe. Se não existir insere ele no grafo
    if (!this->searchNode(id)) {
        this->insertNode(id);
        // cout << "Inserindo " << id << endl;
    }
    // Procura se o nó target_id existe. Se não existir insere ele no grafo
    if (!this->searchNode(target_id)) {
        this->insertNode(target_id);
        // cout << "Inserindo " << target_id << endl;
    }

    Vertice *nodeId = this->getNode(id);
    Vertice *nodeTargetId = this->getNode(target_id);

    if (this->getDirected()) {
        // Cria a aresta => id -> target_id
        nodeId->insertEdge(target_id, weight);

        // Aumenta os graus de saída e de entrada
        nodeId->incrementOutDegree();
        nodeTargetId->incrementInDegree();
    } else {
        // Cria a aresta => id - target_id
        if (!this->searchEdge(id, target_id)) {
            nodeId->insertEdge(target_id, weight);
            nodeTargetId->insertEdge(id, weight);
            // Aumenta os graus de saída e de entrada
            nodeId->incrementOutDegree();
            nodeTargetId->incrementOutDegree();
            nodeId->incrementInDegree();
            nodeTargetId->incrementInDegree();
        }
    }
    this->number_edges++;
}

bool Grafo::searchEdge(int id, int target_id) {
    Vertice *node = this->getNode(id);
    Vertice *targetNode = this->getNode(target_id);
    Aresta *edge = nullptr;

    edge = node->getFirstEdge();
    while (edge != nullptr) {
        if (edge->getTargetId() == target_id) {
            return true;
        }
        edge = edge->getNextEdge();
    }

    return false;
}

void Grafo::removeNode(int id) {}

bool Grafo::searchNode(int id) {
    Vertice *node = this->getFirstNode();
    while (node != nullptr) {
        if (node->getId() == id)
            return true;
        node = node->getNextNode();
    }
    return false;
}

Vertice *Grafo::getNode(int id) {
    Vertice *node = this->getFirstNode();

    while (node != nullptr) {
        if (node->getId() == id)
            return node;
        node = node->getNextNode();
    }
    return nullptr;
}
void Grafo::printGraph() {
    Vertice *node = this->getFirstNode();
    Aresta *edge = nullptr;

    cout << "Lista de adjacência" << endl;

    if (this->getDirected()) {
        if (this->getWeightedEdge()) {
            while (node != nullptr) {
                edge = node->getFirstEdge();
                cout << node->getId() << " -> ";

                while (edge != nullptr) {
                    cout << edge->getTargetId() << " -(" << edge->getWeight() << ")-> ";
                    edge = edge->getNextEdge();
                }

                cout << "null" << endl;
                node = node->getNextNode();
            }
        } else {
            while (node != nullptr) {
                edge = node->getFirstEdge();
                cout << node->getId() << " -> ";

                while (edge != nullptr) {
                    cout << edge->getTargetId() << " -> ";
                    edge = edge->getNextEdge();
                }

                cout << "null" << endl;
                node = node->getNextNode();
            }
        }
    } else {
        if (this->getWeightedEdge()) {
            while (node != nullptr) {
                edge = node->getFirstEdge();
                cout << node->getId() << " - ";

                while (edge != nullptr) {
                    cout << edge->getTargetId() << " -(" << edge->getWeight() << ")- ";
                    edge = edge->getNextEdge();
                }

                cout << "null" << endl;
                node = node->getNextNode();
            }
        } else {
            while (node != nullptr) {
                edge = node->getFirstEdge();
                cout << node->getId() << " -> ";

                while (edge != nullptr) {
                    cout << edge->getTargetId() << " - ";
                    edge = edge->getNextEdge();
                }

                cout << "null" << endl;
                node = node->getNextNode();
            }
        }
    }
}

void Grafo::printGraphDot(ofstream &file) {
    if (file.is_open()) {
        cout << "Salvando o grafo" << endl;

        Vertice *node = this->getFirstNode();
        Aresta *edge;

        vector<bool> visited;

        while (node != nullptr) {
            visited.push_back(false);
            node = node->getNextNode();
        }
        node = this->getFirstNode();
        // Verifica se é ou não direcionado
        if (this->getDirected()) {
            file << "digraph { \n";
        } else {
            file << "strict graph { \n";
        }

        // Verifica se o nó tem peso
        if (this->getWeightedNode() == 1) {
            while (node != nullptr) {
                file << "   " << node->getId() << " [weight = ";
                file << node->getWeight() << "] \n";
                node = node->getNextNode();
            }
        }

        node = this->getFirstNode();

        while (node != nullptr) {
            if(!visited.at(node->getId())) {
                visited.at(node->getId()) = true;
                edge = node->getFirstEdge();
                while (edge != nullptr) {
                    if(!visited.at(edge->getTargetId())) {
                        file << "   " << node->getId();
                        if (this->getDirected()) {
                            file << " -> ";
                        } else {
                            file << " -- ";
                        }
                        file << edge->getTargetId();

                        if (this->getWeightedEdge()) {
                            file << " [label=" << edge->getWeight();
                            file << ",weight=" << edge->getWeight() << "]";
                        }
                        file << "\n";
                    }
                    edge = edge->getNextEdge();
                }
            }
            node = node->getNextNode();
        }

        file << "} \n";
    } else {
        cout << "Falha ao abrir o arquivo";
    }
}

vector<Grafo*> Grafo::guloso(bool random, double *result, float alfa) {
    vector<Grafo*> solution;
    *result = 0;
    vector<bool> visitedNodes;
    int countVisitedNodes = 0;
    vector<vector<bool>> visitedEdges;

    visitedEdges.resize(this->getOrder());
    for(int i = 0; i < this->getOrder(); i++) {
        visitedEdges.at(i).resize(this->getOrder());
        visitedNodes.push_back(false);
    }

    float resultBenefit = 0;

    for(int i = 0; i < this->cluster; i++) {
        tuple<int, int> limits = this->clustersLimits.at(i);
        Grafo *cluster = new Grafo(get<0>(limits), get<1>(limits));
        solution.push_back(cluster);
    }

    for(int i = 0; i < this->cluster; i++) {
        Vertice *node;

        if(random) {
            node = this->verticeValido(0.0f, alfa * (this->getOrder()- 1));
            // position = (int)(rand() % (int)(alfa*this->getOrder()));
        } else {
            node = this->getNode(i);
        }
        int position = node->getId();
        // cout << position << endl;
        

        if(node == nullptr || visitedNodes.at(position) == true) {
            i--;
            continue;
        }

        visitedNodes.at(position) = true;
        countVisitedNodes++;

        // solution.at(i)->insertNodeAndWeight(node->getId(), node->getWeight());
        solution.at(i)->insertNodeAndWeight(position, node->getWeight());
        solution.at(i)->currentLimit += node->getWeight();
    }

    priority_queue<tuple<float, int, int>> candidates;
    
    for(int i = 0; i < this->cluster; i++) {
        Grafo *clusterGraph = solution.at(i);
        int auxId = clusterGraph->getFirstNode()->getId();
        
        while(clusterGraph->currentLimit < clusterGraph->inferiorLimit || candidates.empty()) {
            
            for(int j = 0; j < this->getOrder(); j++) {
                if(!visitedNodes.at(j)) {
                    float auxDistance = distancia(auxId, j);
                    tuple<float, int, int> candidate(auxDistance, auxId, j);
                    candidates.push(candidate);
                }
            }
            
            tuple<float, int, int> candidate = candidates.top();
            float distance = get<0>(candidate);
            tuple<int, int> twoNodes(get<1>(candidate), get<2>(candidate));
            candidates.pop();

            Vertice *graphNode1 = clusterGraph->getNode(get<0>(twoNodes));
            Vertice *graphNode2 = this->getNode(get<1>(twoNodes));

            if(graphNode1 == nullptr) {
                graphNode1 = clusterGraph->getNode(get<1>(twoNodes));
                graphNode2 = this->getNode(get<0>(twoNodes));
            }
            if(
                clusterGraph->currentLimit + graphNode2->getWeight() < clusterGraph->upperLimit && 
                visitedNodes.at(graphNode2->getId()) == false
            ) {
                clusterGraph->insertNodeAndWeight(graphNode2->getId(), graphNode2->getWeight());
                clusterGraph->maxBenefit += distance;
                visitedEdges.at(graphNode2->getId()).at(graphNode1->getId()) = true;
                visitedEdges.at(graphNode1->getId()).at(graphNode2->getId()) = true;

                Vertice *clusterNode = clusterGraph->getFirstNode();
                while(clusterNode != nullptr) {
                    if(
                        visitedEdges.at(graphNode2->getId()).at(graphNode1->getId()) == false && 
                        visitedEdges.at(graphNode1->getId()).at(graphNode2->getId()) == false
                    ) {
                        float auxDistance = distancia(graphNode2->getId(), clusterNode->getId());
                        clusterGraph->maxBenefit += auxDistance;
                        resultBenefit += auxDistance;
                    }
                    clusterNode = clusterNode->getNextNode();
                }
                resultBenefit += distance;

                clusterGraph->currentLimit += graphNode2->getWeight();
                visitedNodes.at(graphNode2->getId()) = true;
                countVisitedNodes++;
                auxId = graphNode2->getId();
            }
        }
    }

    for(int i = 0; i < this->getOrder(); i++) {
        for(int j = 0; j < this->getOrder(); j++) {
            if(visitedNodes.at(i) == false || visitedNodes.at(j) == false) {
                float auxDistance = distancia(i, j);
                tuple<float, int, int> candidate(auxDistance, i, j);
                candidates.push(candidate);
            }
        }
    }

    while(countVisitedNodes < this->getOrder() && !candidates.empty()) {
        tuple<float, int, int> candidate = candidates.top();
        float distance = get<0>(candidate);
        tuple<int, int> twoNodes(get<1>(candidate), get<2>(candidate));
        candidates.pop();

        if(
            !(visitedNodes.at(get<0>(twoNodes)) == true && visitedNodes.at(get<1>(twoNodes)) == true) && 
            !(visitedNodes.at(get<0>(twoNodes)) == false && visitedNodes.at(get<1>(twoNodes)) == false)
        ) {
            for(int i = 0; i < this->cluster; i++) {
                Grafo *cluster = solution.at(i);

                Vertice *graphNode1 = cluster->getNode(get<0>(twoNodes));
                Vertice *graphNode2 = this->getNode(get<1>(twoNodes));

                if(graphNode1 == nullptr) {
                    graphNode1 = cluster->getNode(get<1>(twoNodes));
                    graphNode2 = this->getNode(get<0>(twoNodes));
                }

                if(
                    cluster->currentLimit + graphNode2->getWeight() <= cluster->upperLimit &&
                    visitedNodes.at(graphNode2->getId()) == false
                ) {
                    cluster->insertNodeAndWeight(graphNode2->getId(), graphNode2->getWeight());
                    cluster->maxBenefit += distance;
                    resultBenefit += distance;

                    visitedEdges.at(get<0>(twoNodes)).at(get<1>(twoNodes)) = true;
                    visitedEdges.at(get<1>(twoNodes)).at(get<0>(twoNodes)) = true;

                    Vertice *clusterNode = cluster->getFirstNode();
                    while(clusterNode != nullptr) {
                        if(
                            visitedEdges.at(graphNode2->getId()).at(clusterNode->getId()) == false && 
                            visitedEdges.at(clusterNode->getId()).at(graphNode2->getId()) == false
                        ) {
                            float auxDistance = distancia(graphNode2->getId(), clusterNode->getId());
                            cluster->maxBenefit += auxDistance;
                            resultBenefit += auxDistance;
                        }
                        clusterNode = clusterNode->getNextNode();
                    }
                    cluster->currentLimit += graphNode2->getId();
                    visitedNodes.at(graphNode2->getId()) = true;
                    countVisitedNodes++;
                }
            }
        }
    }

    *result = resultBenefit;
    return solution;
}

void Grafo::agmGuloso() {
    auto start = chrono::steady_clock::now();

    double result = 0;
    vector<Grafo*> sol = guloso(0, &result, 0);

    auto end = chrono::steady_clock::now();
    cout << "Demorou  "
            << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            << " ms para ler o arquivo de entrada." << endl;
    // cout << "Qualidade Solucao: " << qualidadeSolucao(result) << "%" << endl;
    float resLitera = 225003.70;
    cout << "Qualidade Obtida: " << result << " | Qualidade Literatura: " <<  resLitera << endl;
    cout << "Diferença (%): " << (result/resLitera) << endl;
    cout << "Diferença (%): " << (resLitera/result) << endl;
    if (result > 0) {
        cout << "Conseguiu alguma solucao viavel" << endl;
    } else {
        cout << "Nao conseguiu nenhuma solucao viavel" << endl;
    }

    imprimeCluster(sol, 2, result);

    // imprimeCluster(sol, 2, result);
    // output("AlgoritmoGuloso.txt", sol, qualidadeSolucao(result));
}
// GULOSOS


void Grafo::gulosoRandAdap(){
    auto start = chrono::steady_clock::now();

    float melhor = 0;
    double resultado = 0;
    int criterio_parada=100;
    float cof_randomizacao;

    // cout << "Escolha um coeficiente de randomizacao: " << endl;
    // cin >> cof_randomizacao;
    for (int i = 0; i < rand(); i++)
        cof_randomizacao = 0 + (float) (rand()) / ((float) (RAND_MAX / (1 - 0)));

    vector<Grafo *> solution, best_solution;

    cout << "Coeficiente de randomizacao: " << cof_randomizacao << endl;

    int i=0;
    while(i < criterio_parada) {
        solution = guloso(1, &resultado, cof_randomizacao);
        // cout << "guloso concluido" << endl;
        if (resultado > melhor) {
            melhor = resultado;
            best_solution =solution;
        }
        i++;
    }
    cout << std::setprecision(2) << std::fixed;
    auto end = chrono::steady_clock::now();

    cout << "Beneficio da melhor solucao: " << melhor <<endl;
    if (melhor > 0) {
        cout << "O guloso randomizado obteve alguma solucao viavel" << endl;
    } else {
        cout << "O guloso randomizado reativo nao obteve nenhuma solucao viavel" << endl;
    }

    //output("AlgoritmoGulosoRandomizadoAdaptativo.txt", melhorSol, qualidadeSolucao(maior));
}
struct media {
    float soma;
    float numSolucoes;
    float media;
};
void inicializaVetores(vector<float> &prob, vector<media> &medias, size_t numAlfas) {
    media aux{0, 0, 1};
    auto auxNumAlfas = numAlfas;
    float auxProb = 1.0f / static_cast<float>(auxNumAlfas);
    for (int i = 0; i < numAlfas; i++) {
        prob.push_back(auxProb);
        medias.push_back(aux);
    }
}
void atualizaProbabilidades(vector<media> &medias, vector<float> &prob, vector<float> &solBest, vector<float> &q) {
    float somaQ = 0;
    float melhorSolucao = *(max_element(solBest.begin(), solBest.end()));
    for (int i = 0; i < medias.size(); i++) {
        q[i] = pow((melhorSolucao / medias[i].media), 2);
        somaQ += q[i];
    }
    for (int i = 0; i < medias.size(); i++) {
        prob[i] = q[i] / somaQ;
    }
}
float escolheAlfa(vector<float> &prob, vector<float> &alfas) {
    float soma = 0;
    int r = rand() % 101;
    for (int i = 0; i < prob.size(); i++) {
        soma += (prob[i] * 100);
        if (soma >= r) {
            return alfas[i];
        }
    }
    return alfas[alfas.size() - 1];
}
void atualizaMedias(vector<media> &medias, float solucao, vector<float> &alfas, float alfa) {
    size_t aux = 0;
    float auxSolucao;
    auxSolucao = solucao;
    for (size_t i = 0; i < alfas.size(); i++) {
        if (alfa == alfas[i]) {
            aux = i;
            break;
        }
    }
    medias[aux].soma = medias[aux].soma + auxSolucao;
    medias[aux].numSolucoes++;
    medias[aux].media = medias[aux].soma / medias[aux].numSolucoes;
}


void Grafo::gulosoReativo() {
    auto start = chrono::steady_clock::now();

    vector<float> alfas{0.05f, 0.10f, 0.15f, 0.30f, 0.50f}, solBest, probabilidade, q;
    vector<media> medias;
    vector<Grafo *> solution, melhorSol;
    int criterio_parada = 2500;
    int numBloco = 50;
    float solucao, bestBenefit= 0;


    for (int i = 0; i < alfas.size(); i++) {
        q.push_back(0.00f);
        solBest.push_back(0.00f);
    }
    inicializaVetores(probabilidade, medias, alfas.size());
    for (int i = 1; i <= criterio_parada; i++) {
        if (i % numBloco == 0) {
            atualizaProbabilidades(medias, probabilidade, solBest, q);
        }
        float cof_randomizado = escolheAlfa(probabilidade, alfas);
        /**/
        double result;
        float semente;
        cout << "Escolha um coeficiente de randomizacao: " << cof_randomizado << endl;
        solution = guloso(1, &result, cof_randomizado);
        cout << "Beneficio: " << result << endl
             << endl;
        if (result > bestBenefit) {
            bestBenefit = result;
            melhorSol = solution;
            solBest[cof_randomizado] = bestBenefit;
        }
        for (auto &i : solution) {
            delete i;
        }
        atualizaMedias(medias, bestBenefit, alfas, cof_randomizado);
    }
    float auxSolBest = 0;
    for (int i = 0; i < solBest.size(); i++) {
        if (solBest[i] > auxSolBest) {
            auxSolBest = solBest[i];
        }
    }
    cout << std::setprecision(2) << std::fixed;
    auto end = chrono::steady_clock::now();
    cout << "Beneficio da Melhor Solucao: " << auxSolBest << endl;
    // cout << "Semente da melhor solucao: " << alfas[auxSolBest] << endl;
    // cout << "Qualidade Solucao: " << qualidadeSolucao(auxSolBest) << "%" << endl;
    if (auxSolBest > 0) {
        cout << "Conseguiu alguma solucao viavel" << endl;
    } else {
        cout << "Nao conseguiu nenhuma solucao viavel" << endl;
    }
    //output("AlgoritmoGulosoRandomizadoReativo.txt", melhorSol, qualidadeSolucao(auxSolBest));
}

float Grafo::distancia(int node1, int node2) {
    // START - Get Distance between two nodes
    Vertice *nodeAux1 = this->getNode(node1);
    Vertice *nodeAux2 = this->getNode(node2);
    Aresta *edgeAux;
    float auxDistance = 0;
    if(nodeAux1 != nullptr) {
        edgeAux = nodeAux1->getFirstEdge();
        while (edgeAux != nullptr) {
            if(edgeAux->getTargetId() == nodeAux2->getId())
                auxDistance = edgeAux->getWeight();
            edgeAux = edgeAux->getNextEdge();
        }
    }
    // END - Get Distance between two nodes
    return auxDistance;
}


void Grafo::qualidade(float result) {
    switch (this->fileType) {
        case 1: {
            // RanReal240_01.txt
            float literatureResult = 225003.70;
            break;
        }
        case 2: {
            // RanReal240_04.txt
            float literatureResult = 225683.17;
            break;
        }
        case 3: {
            // RanReal240_07.txt
            float literatureResult = 209305.70;
            break;
        }
        case 4: {
            // RanReal480_01.txt
            float literatureResult = 556126.86;
            break;
        }
        case 5: {
            // RanReal480_04.txt
            float literatureResult = 522790.22;
            break;
        }
        case 6: {
            // RanReal960_01.30.txt
            float literatureResult = 1340369.47;
            break;
        }
        case 7: {
            // Sparse82_02.txt
            float literatureResult = 1306.64;
            break;
        }
        case 8: {
            // 20_5_270001
            float literatureResult = 0;
            break;
        }
        case 9: {
            // 20_10_270001
            float literatureResult = 2148.00;
            break;
        }
        case 10: {
            // 30_5_270003
            float literatureResult = 920.00;
            break;
        }
        default: {
            cout << "Exit!!!" << endl;
        }

    }
}



Vertice* Grafo::verticeValido(float min, float max) {
    float random = ((float)rand()) / (float)RAND_MAX;
    float diff = max - min;
    float r = random * diff;
    int idRandom = min + r;
    // find the node with at idRandom
    list<int> candidates;
    for (int i = 0; i < this->getOrder(); i++) {
        candidates.push_back(i);
    }

    int i = 0;
    for (auto it = candidates.begin(); it != candidates.end(); ++it)
    {
        if (i == idRandom)
        {
            auto node = this->getNode(*it);
            candidates.erase(it);
            return node;
        }
        i++;
    }

    return nullptr;
}


void Grafo::imprimeCluster(vector<Grafo *> solucao, int option, float resultBeneficio)
{
    float totalBeneficio = 0.0;

    for (int i = 0; i < this->cluster; i++)
    {
        Grafo *cluster = solucao[i];

        if (option == 2)
        {
            cout << "----------IMPRIME CLUSTER " << i + 1 << "----------" << endl;
            cout << "Beneficio " << cluster->maxBenefit << endl;
            totalBeneficio += cluster->maxBenefit;
        }

        if (option == 1)
        {
            cluster->imprimirVertices();
        }
        else if (option == 2)
        {
            cluster->imprimirVertices2();
        }

        cout << endl;
    }

    if (option == 2)
    {
        cout << std::setprecision(2) << std::fixed;
    }
    cout << "\n\nBeneficio final: " << totalBeneficio << endl;
}

void Grafo::imprimirVertices2()
{
    Vertice *node = this->first_node;
    int cont = 0;

    cout << "Limite | " << this->inferiorLimit << " <= " << this->currentLimit << " <= " << this->upperLimit << ""
         << endl;

    while (node != nullptr)
    {
        cout << node->getId() << ",";
        node = node->getNextNode();
        cont++;
    }
}

void Grafo::imprimirVertices()
{
    Vertice *node = this->first_node;
    int cont = 0;

    while (node != nullptr)
    {
        cout << node->getId() << ",";
        node = node->getNextNode();
        cont++;
    }
}

void Grafo::output(string output_file, vector<Grafo*> solution, float quality){
    fstream output_File(output_file, ios::in);

    output_File << "Clusterização: " << endl;

    for (int i = 0; i < this->cluster; i++) {
        Grafo *cluster = solution[i];
        output_File << "---------- CLUSTER ----------" << endl;

        output_File << "Beneficio do Cluster: " << cluster->maxBenefit << endl;

        Vertice *node = cluster->getFirstNode();

        output_File << "Limite: " << cluster->inferiorLimit << " <->"/*<< cluster->getLimit << "<->"*/ << cluster->upperLimit << endl;

        output_File << "Nós do Cluster: " << endl;
        while (node != nullptr) {
            output_File << node->getId() << ";";
            node = node->getNextNode();
        }

        output_File << endl;
    }
    output_File << "Qualidade: " << quality << "%" << endl;

    cout << "O arquivo " << output_file << " foi gerado com sucesso.";
}

