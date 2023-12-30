#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define MAX_CMD 19
#define AGG_STAZ "aggiungi-s"   // aggiungi-stazione
#define DEM_STAZ "d"            // demolisci-stazione
#define AGG_AUTO "aggiungi-a"   // aggiungi-auto
#define ROT_AUTO "r"            // rottama-auto
#define PIA_PERC "p"            // pianifica-percorso
#define AGG_STAZ_DIM 10
#define DEM_STAZ_DIM 1
#define AGG_AUTO_DIM 10
#define ROT_AUTO_DIM 1
#define PIA_PERC_DIM 1

#define NIL NULL

#define AGGIUNTA "aggiunta"
#define NON_AGGI "non aggiunta"
#define DEMOLITA "demolita"
#define NON_DEMO "non demolita"
#define ROTTAMAT "rottamata"
#define NON_ROTT "non rottamata"
#define NESS_PER "nessun percorso"
#define AGGIUNTA_DIM 8
#define NON_AGGI_DIM 12
#define DEMOLITA_DIM 8
#define NON_DEMO_DIM 12
#define ROTTAMAT_DIM 9
#define NON_ROTT_DIM 13
#define NESS_PER_DIM 15

#define N_AGGIUNGI_AUTO_PAR 2

#define MAX_INT 2147483647

#define CHECK 0
#define DEBUG 0

typedef struct {

    struct leaf* left;
    struct leaf* right;
    struct leaf* p;
    int key;
    int maxAutonomy;
    struct leaf* veichles;
    int distance;
    struct leaf* prev;
    struct leaf* next;

}leaf;

int readCommand(char**);
int stringCompare(char*, char*, int);
void readAggiungiStazioneParameters(int**);
void insertInBST(leaf**, leaf*);
leaf* searchInBST(leaf*, int);
void printStdoutOptimized(char*, int);
int hashFunction(int);
int readDemolisciStazioneParameter();
void freeBST(leaf*);
void removeFromBST(leaf**, leaf*);
leaf* nextInBST(leaf*);
leaf* minInBST(leaf*);
void readTwoIntegerParameters(int**);
void addVeichle(leaf*, leaf*);
int scrapVeichle(leaf*, int);
leaf* maxInBST(leaf*);
leaf* previousInBST(leaf*);
void printPianificaPercorso(int*, int);
void printIntStdoutOptimized(int);
void enqueueWithPrio(leaf**, leaf*);
leaf* removeMinFromQueue(leaf**);
void dijkstraForward(int, leaf*, int, leaf*, int**, int*, int*, int);

void inOrderBST(leaf*);
void inOrderBST_STAR(leaf*);

int main() {

    int end = 0;
    char* command = NULL;
    int* commandArguments = NULL;
    leaf* distanceBst = NIL;
    int* output = NULL;
    int num_stations = 0;
    int max_station = 0;

    // start reading from stdin
    do{
        if(!readCommand(&command)) {    // check if the input file is empty
            free(command);  // free all the allocated pointers
            free(commandArguments);
            freeBST(distanceBst);
            end = 1;    // set flag to close software
        }
        else {
            if(stringCompare(command, AGG_STAZ, AGG_STAZ_DIM)) {    // check what type of command
                readAggiungiStazioneParameters(&commandArguments);  // read arguments for aggiungi-stazione
                if(commandArguments[0] > max_station) {
                    max_station = commandArguments[0];
                }
                leaf* found = searchInBST(distanceBst, commandArguments[0]);    // check if already exist a station
                if(found == NIL) {
                    leaf* distanceLeaf = (leaf*) malloc(sizeof(leaf));  // if not found add it to the stations BST
                    distanceLeaf->key = commandArguments[0];
                    distanceLeaf->left = NIL;
                    distanceLeaf->right = NIL;
                    distanceLeaf->maxAutonomy = 0;
                    distanceLeaf->veichles = NIL;
                    if(commandArguments[1] > 0) {
                        for(int i = 2; i < commandArguments[1] + 2; i++) {
                            leaf* veichleLeaf = (leaf*) malloc(sizeof(leaf));
                            veichleLeaf->key = commandArguments[i];
                            veichleLeaf->left = NIL;
                            veichleLeaf->right = NIL;
                            veichleLeaf->p = NIL;
                            veichleLeaf->veichles = NIL;
                            veichleLeaf->maxAutonomy = 0;
                            insertInBST((leaf**) &(distanceLeaf->veichles), veichleLeaf);
                            if(commandArguments[i] > distanceLeaf->maxAutonomy) {
                                distanceLeaf->maxAutonomy = commandArguments[i];
                            }
                        }
                    }
                    insertInBST(&distanceBst, distanceLeaf);
                    printStdoutOptimized(AGGIUNTA, AGGIUNTA_DIM);    // print aggiunta
                    num_stations++;
                }
                else {
                    printStdoutOptimized(NON_AGGI, NON_AGGI_DIM);    // if already exist, print non aggiunta
                }
            }
            else if(stringCompare(command, DEM_STAZ, DEM_STAZ_DIM)) {
                int toDemolish = readDemolisciStazioneParameter();
                leaf* found = searchInBST(distanceBst, toDemolish);
                if(found != NULL) {
                    removeFromBST(&distanceBst, found);
                    printStdoutOptimized(DEMOLITA, DEMOLITA_DIM);
                    num_stations--;
                }
                else {
                    printStdoutOptimized(NON_DEMO, NON_DEMO_DIM);
                }
            }
            else if(stringCompare(command, AGG_AUTO, AGG_AUTO_DIM)) {
                readTwoIntegerParameters(&commandArguments);
                leaf* found = searchInBST(distanceBst, commandArguments[0]);
                if(found != NULL) {
                    leaf* veichleToAdd = malloc(sizeof(leaf));
                    veichleToAdd->key = commandArguments[1];
                    veichleToAdd->p = NIL;
                    veichleToAdd->left = NIL;
                    veichleToAdd->right = NIL;
                    veichleToAdd->maxAutonomy = 0;
                    veichleToAdd->veichles = NIL;
                    addVeichle(found, veichleToAdd);
                    printStdoutOptimized(AGGIUNTA, AGGIUNTA_DIM);
                }
                else {
                    printStdoutOptimized(NON_AGGI, NON_AGGI_DIM);
                }
            }
            else if(stringCompare(command, ROT_AUTO, ROT_AUTO_DIM)) {
                readTwoIntegerParameters(&commandArguments);
                leaf* found = searchInBST(distanceBst, commandArguments[0]);
                if(found != NULL) {
                    if(scrapVeichle(found, commandArguments[1])) {
                        printStdoutOptimized(ROTTAMAT, ROTTAMAT_DIM);
                    }
                    else {
                        printStdoutOptimized(NON_ROTT, NON_ROTT_DIM);
                    }
                }
                else {
                    printStdoutOptimized(NON_ROTT, NON_ROTT_DIM);
                }
            }
            else if(stringCompare(command, PIA_PERC, PIA_PERC_DIM)) {
                readTwoIntegerParameters(&commandArguments);
                if(commandArguments[0] == commandArguments[1]) {
                    printIntStdoutOptimized(commandArguments[0]);
                    putchar_unlocked('\n');
                }
                else if(commandArguments[0] < commandArguments[1]) {
                    leaf* startLeaf = searchInBST(distanceBst, commandArguments[0]);
                    output = (int*) malloc(sizeof(int) * num_stations);
                    int exist = 0;
                    int count = 0;
                    dijkstraForward(num_stations, startLeaf, commandArguments[1], distanceBst, &output, &count, &exist, max_station+1);
                    if(!exist) {
                        printStdoutOptimized(NESS_PER, NESS_PER_DIM);
                    }
                    else {
                        printPianificaPercorso(output, count);
                    }
                    // free the output vector mallocated
                    free(output);
                }
                else {
                    printf("PIANIFICA PERCORSO BACKWARDS\n");
                    /*leaf* startLeaf = searchInBST(distanceBst, commandArguments[0]);
                    output = (int*) malloc(sizeof(int) * num_stations);
                    int exist = 0;
                    int count = 0;
                    dijkstraBackwards(num_stations, startLeaf, commandArguments[1], hashtable, &output, &count, &exist, max_station+1);
                    if(!exist) {
                        printStdoutOptimized(NESS_PER, NESS_PER_DIM);
                    }
                    else {
                        printPianificaPercorso(output, count);
                    }
                    // free the output vector mallocated
                    free(output);*/
                }
            }
        }
    }while(!end);

    return 0;
}

int readCommand(char** cmd) {   // read command from stdin
    char read = '\0';
    int count = 0;

    do {
        read = getchar_unlocked();  // read char from stdin
        if(read == EOF) {   // chek for end of file
            return 0;
        }
        else if(read != ' ') {
            if(count == 0) {    // if first char readed allocate memory
                if(*cmd != NULL) {
                    free(*cmd);
                }
                (*cmd) = (char*) malloc(sizeof(char) * MAX_CMD);
            }
            (*cmd)[count] = read;   // save readed char
            count++;
        }
    } while(read != ' ');   // stop when read a space
    return 1;
}

int stringCompare(char* a, char* b, int dim) {  // compare two strings
    for(int i = 0; i < dim; i++) {
        if(a[i] != b[i]) {  // if chars are not equals return false
            return 0;
        }
    }
    return 1;   // return true
}

void readAggiungiStazioneParameters(int** vector) { // read all the parameters of the command "aggiungi-stazione"
    
    char read = '\0';
    int distance = 0, size = 0, val;

    
    do{ // read from stdin distance
        read = getchar_unlocked();
        if(read != ' ') {
            distance = distance*10 + read - '0';
        }
    }while(read != ' ');
    do{ // read from stdin # of veichles
        read = getchar_unlocked();
        if(read != ' ' && read != '\n' && read != EOF) {
            size = size*10 + read - '0';
        }
    }while(read != ' ' && read != '\n' && read != EOF);
    if(*vector != NULL) {
        free(*vector);
    }
    (*vector) = (int*) malloc(sizeof(int) * (size + 2));    // allocate memory
    (*vector)[0] = distance;    // save as first value # veichles
    (*vector)[1] = size;
    if(size > 0) {
        for(int i = 2; i < size + 2; i++) { // read all veichle values from stdin
            read = '\0';
            val = 0;
            do {
                read = getchar_unlocked();
                if(read != ' ' && read != '\n' && read != EOF) {
                    val = val*10 + read - '0';
                }
            }while(read != ' ' && read != '\n' && read != EOF);
            (*vector)[i] = val;
        }
    }
    return;
}

void insertInBST(leaf** T, leaf* x) {   // insert value into a BST

    leaf* pre = NIL;
    leaf* cur = (*T);

    while(cur != NIL) {
        pre = cur;
        if(x->key < cur->key) {
            cur = (leaf*) cur->left;
        }
        else {
            cur = (leaf*) cur->right;
        }
    }
    x->p = (struct leaf*) pre;
    if(pre == NIL) {
        (*T) = x;
    }
    else if(x->key < pre->key) {
        pre->left = (struct leaf*) x;
    }
    else {
        pre->right = (struct leaf*) x;
    }
    
    return;
}

leaf* searchInBST(leaf* T, int x) { // search for value in BST

    if(T == NIL || (T != NIL && T->key == x)) {
        return T;
    }
    if(T->key < x) {
        return searchInBST((leaf*) T->right, x);
    }
    else {
        return searchInBST((leaf*) T->left, x);
    }

}

void printStdoutOptimized(char* str, int dim) {
    
    for (int i = 0; i < dim; i++) {
        putchar_unlocked(str[i]);
    } 
    putchar_unlocked('\n');
    return;
    
}

int abs(int x) {

    return x < 0 ? (-1) * x : x;

}

int readDemolisciStazioneParameter() {  // read distance parameter of the station to demolish

    char read = '\0';
    int val = 0;

    do {    
        read = getchar_unlocked();  
        if(read != ' ' && read != '\n' && read != EOF) {
            val = val*10 + read - '0';
        }
    }while(read != ' ' && read != '\n' && read != EOF);
    return val;

}

void freeBST(leaf* T) { // post order BST visit to free all the leafs
    
    if(T->left != NIL) {
        freeBST((leaf*) T->left);
    }
    if(T->right != NIL) {
        freeBST((leaf*) T->right);
    }
    if(T->veichles != NIL) {
        freeBST((leaf*) T->veichles);
    }
    free(T);
    return;

}

void removeFromBST(leaf** T, leaf* x) { // remove element from BST

    leaf* toDelete = NIL;
    leaf* under = NIL;

    if(x->left == NIL || x->right == NIL) {
        toDelete = x;
    }
    else {
        toDelete = nextInBST(x);
    }
    if(toDelete->left != NIL) {
        under = (leaf*) toDelete->left;
    }
    else {
        under = (leaf*) toDelete->right;
    }
    if(under != NIL) {
        under->p = toDelete->p;
    }
    leaf* toDelete_p = (leaf*) toDelete->p;
    if(toDelete->p == NIL) {
        (*T) = under;
    }
    else if(toDelete == (leaf*) toDelete_p->left) {
        toDelete_p->left = (struct leaf*) under;
    }
    else {
        toDelete_p->right = (struct leaf*) under;
    }
    if(toDelete != x) {
        x->key = toDelete->key;
        x->maxAutonomy = toDelete->maxAutonomy;
        if(x->veichles != NIL) {
            freeBST((leaf*) x->veichles);
        }
        x->veichles = toDelete->veichles;
    }
    else {
        if(x->veichles != NIL) {
            freeBST((leaf*) x->veichles);
        }
    }
    free(toDelete);
    return;

}

leaf* minInBST(leaf* T) {   // search for element with min key in BST

    leaf* curr = T;

    while(curr->left != NIL) {
        curr = (leaf*) curr->left;
    }
    return curr;

}

leaf* nextInBST(leaf* x) {  // search element with next value of the key in BST

    leaf* y = NIL;

    if(x->right != NIL) {
        return minInBST((leaf*) x->right);
    }
    y = (leaf*) x->p;
    while(y != NIL && (leaf*) y->right == x) {
        x = y;
        y = (leaf*) y->p;
    }
    return y;

}

void readTwoIntegerParameters(int** vector) {  // read arguments for aggiungi-auto cmd

    char read = '\0';
    int val = 0;

    if(*vector != NULL) {
        free(*vector);
    }
    (*vector) = (int*) malloc(sizeof(int) * N_AGGIUNGI_AUTO_PAR);
    for(int i = 0; i < N_AGGIUNGI_AUTO_PAR; i++) {
        val = 0;
        do {
            read = getchar_unlocked();  // read two integer values
            if(read != ' ' && read != '\n' && read != EOF) {
                val = val*10 + read - '0';
            }
        }while(read != ' ' && read != '\n' && read != EOF);
        (*vector)[i] = val;
    }
    return;

}

void addVeichle(leaf* T, leaf* veichleToAdd) {   // adds a new veichle into the BST of veichles of one station

    insertInBST((leaf**) &(T->veichles), veichleToAdd);   // add new veichle in BST
    if(veichleToAdd->key > T->maxAutonomy) {
        T->maxAutonomy = veichleToAdd->key;
    }

    return;

}

int scrapVeichle(leaf* T, int autonomy) {    // remove, if exists, an element from the BST

    leaf* toScrap = searchInBST((leaf*) T->veichles, autonomy);  // check if exists the element to scrap
    if(toScrap == NULL) {   // if not exists return false
        return 0;
    }
    else {
        int toScrapAutonomy = toScrap->key;         // save the distance
        removeFromBST((leaf**) &(T->veichles), toScrap);    // else remove it and return true
        if(T->veichles == NULL) {
            T->maxAutonomy = 0;
        }
        else if(T->maxAutonomy == toScrapAutonomy) {
            T->maxAutonomy = maxInBST((leaf*) T->veichles)->key;
        }
        return 1;
    }

}

leaf* maxInBST(leaf* T) {

    leaf* curr = T;

    while(curr->right != NIL) {
        curr = (leaf*) curr->right;
    }
    return curr;

}

leaf* previousInBST(leaf* x) {  // search element with next value of the key in BST

    leaf* y = NIL;

    if(x->left != NIL) {
        return maxInBST((leaf*) x->left);
    }
    y = (leaf*) x->p;
    while(y != NIL && (leaf*) y->left == x) {
        x = y;
        y = (leaf*) y->p;
    }
    return y;

}

void printIntStdoutOptimized(int x) {

    int N = x, rev, count = 0;
    rev = N;
    if (N == 0) {
        putchar_unlocked('0');
        return ;
    }
    while((rev % 10) == 0) { 
        count++;
        rev /= 10;
    }
    rev = 0;
    while(N != 0) {
        rev = (rev<<3) + (rev<<1) + N % 10;
        N /= 10;
    }
    while(rev != 0) {
        putchar_unlocked(rev % 10 + '0');
        rev /= 10;
    }
    while(count--) {
        putchar_unlocked('0');
    }

}

void printPianificaPercorso(int* stations, int count) {

    for(int i = count - 1; i >= 0; i--) {
        printIntStdoutOptimized(stations[i]);
        if(i != 0) {
            putchar_unlocked(' ');
        }
        else {
            putchar_unlocked('\n');
        }
    }
    return;

}

void dijkstraForward(int nStations, leaf* s, int dest, leaf* distanceBst, int** out, int* count, int* exist, int maxDim) {

    leaf* Q = NULL;
    s->distance = 0;
    s->prev = NULL;
    s->next = NULL;
    leaf* v = s;
    while(v->key <= dest) {
        if(v != s) {
            v->next = NULL;
            v->distance = MAX_INT;
            v->prev = NULL;
        }
        enqueueWithPrio(&Q, v);
        if(v->key != dest) {
            v = nextInBST(v);
        }
        else break;
    }
    while(Q != NULL) {
        leaf* u = removeMinFromQueue(&Q);
        v = nextInBST(u);
        if(v->key > u->key + u->maxAutonomy) {
            return;
        }
        while(v->key <= u->key + u->maxAutonomy && v->key <= dest) {
            int ndist = u->distance + 1;
            if(v != NULL && v->distance > ndist ) {
                v->distance = ndist;
                v->prev = (struct leaf*) u;
            }
            if(v->key != dest) {
                v = nextInBST(v);
            }
            else {
                // SONO ARRIVATO A DESTINAZIONE
                *exist = 1;
                leaf* x = v;
                do {
                    (*out)[*count] = x->key;
                    x = (leaf*) x->prev;
                    *count += 1;
                }while(x != NULL);
                return;
            }
        }
    }
    return;

}

void enqueueWithPrio(leaf** Q, leaf* v)  {

    if(*Q == NULL) {
        *Q = v;
    }
    else {
        leaf* curr = *Q;
        leaf* prev = NULL;
        while(curr->next != NULL) {
            if(curr->distance >= v->distance) {
                prev->next = (struct leaf*) v;
                v->next = (struct leaf*) curr;
                return;
            }
            prev = curr;
            curr = (leaf*) curr->next;
        }
        curr->next = (struct leaf*) v; 
    }
    return;

}

leaf* removeMinFromQueue(leaf** Q) {

    leaf* curr = *Q;
    leaf* pred = NULL;
    leaf* min = NULL;
    leaf* predMin = NULL;

    if((*Q)->next == NULL) {
        min = *Q;
        *Q = NULL;
        return min;
    }

    while(curr != NULL) {

        if(min == NULL) {
            min = curr;
        }
        else if(curr->distance < min->distance) {
            min = curr;
            predMin = pred;
        }
        else if(curr->distance == min->distance && curr->key < min->key) {
            min = curr;
            predMin = pred;
        }
        pred = curr;
        curr = (leaf*) curr->next;
    }
    if(predMin == NULL) {
        *Q = (leaf*) min->next;
    }
    else {
        predMin->next = min->next;
    }
    return min;

}

/*void dijkstraBackwards(int nStations, leaf* s, int dest, bucket** hashTable, int** out, int* count, int* exist, int maxDim) {

    queueElement* Q = NULL;
    leaf* v = s;
    queueElement** leavesToQueueElements = malloc(sizeof(queueElement*) * maxDim);
    for(int i = 0; i < maxDim; i++) {
        leavesToQueueElements[i] = NULL;
    }
    while(v->key >= dest) {
        leavesToQueueElements[v->key] = (queueElement*) malloc(sizeof(queueElement));
        leavesToQueueElements[v->key]->element = v;
        leavesToQueueElements[v->key]->next = NULL;
        if(v == s) {
            leavesToQueueElements[v->key]->distance = 0;
            leavesToQueueElements[v->key]->previous = NULL;
        } 
        else {
            leavesToQueueElements[v->key]->distance = MAX_INT;
            leavesToQueueElements[v->key]->previous = NULL;
        }
        enqueueWithPrio(&Q, leavesToQueueElements[v->key]);
        if(v->key != dest) {
            v = previousInBST(v);
        }
        else {
            break;
        }
    }
    while(Q != NULL) {
        queueElement* u = removeMinFromQueue(&Q);
        int uMaxAutonomy = getMaxAutonomyFromDistance(u->element->key, hashTable, hashFunction(u->element->key));
        v = previousInBST(u->element);
        if(v->key < u->element->key - uMaxAutonomy) {
            freeQueueVector(leavesToQueueElements, maxDim);
            return;
        }
        while(v->key >= u->element->key - uMaxAutonomy && v->key >= dest) {
            int ndist = u->distance + 1;
            queueElement* vQEl = leavesToQueueElements[v->key];
            if(vQEl != NULL && vQEl->distance > ndist ) {
                vQEl->distance = ndist;
                vQEl->previous = (struct queueElement*) u;
            }
            if(v->key != dest) {
                v = previousInBST(v);
            }
            else {
                // SONO ARRIVATO A DESTINAZIONE
                *exist = 1;
                queueElement* x = vQEl;
                do {
                    (*out)[*count] = x->element->key;
                    x = (queueElement*) x->previous;
                    *count += 1;
                }while(x != NULL);
                freeQueueVector(leavesToQueueElements, maxDim);
                return;
            }
        }
    }

}*/

void inOrderBST(leaf* T) {

    if(T == NIL) return;
    inOrderBST((leaf*) T->left);
    printf("%d\n", T->key);
    inOrderBST((leaf*) T->right);

    return;

}

void inOrderBST_STAR(leaf* T) {

    if(T == NIL) return;
    inOrderBST_STAR((leaf*) T->left);
    if(T->key == 7733) {
        printf("VEICHLES OF %d: ", T->key);
        inOrderBST((leaf*) T->veichles);
        printf("\n\n");
    }
    inOrderBST_STAR((leaf*) T->right);

    return;

}