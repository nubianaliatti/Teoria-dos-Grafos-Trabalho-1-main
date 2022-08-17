/*
    CABEÇALHO
*/

#ifndef ARESTA_H
#define ARESTA_H 1

class Aresta {
    private:
    int origem;
    int target_id;
    Aresta *next_edge;
    float weight;

public:
    // Constructor
    Aresta(int target_id);
    // Destructor
    ~Aresta();
    // Getters
    int getTargetId();
    int getOrigem();
    Aresta *getNextEdge();
    float getWeight();
    // Setters
    void setNextEdge(Aresta *edge);
    void setWeight(float weight);
    void setOrigem(int origem);
};

#endif // ARESTA_H