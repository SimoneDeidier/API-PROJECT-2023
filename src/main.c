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

#define HASH_DIM 1024 * 16
#define A 0.6180339887498948482045868343656381177203091798057628621354486227

#define N_AGGIUNGI_AUTO_PAR 2

#define MAX_INT 2147483647

#define DEBUG 0
#define MALLOC 0
#define STAMPA_STRUTTURE 0
#define PPTEST 0
#define CHECK_PREV 0
#define TEST 0
#define DA 1

typedef struct {

    struct leaf* left;
    struct leaf* right;
    struct leaf* p;
    int key;

}leaf;

typedef struct {

    struct bucket* next;
    int distance;
    leaf* veichles;
    int maxAutonomy;

}bucket;

typedef struct {

    leaf* element;
    int distance;
    struct queueElement* previous;
    struct queueElement* next;

}queueElement;

int readCommand(char**);
int stringCompare(char*, char*, int);
void readAggiungiStazioneParameters(int**);
void insertInBST(leaf**, leaf*);
leaf* searchInBST(leaf*, int);
void printStdoutOptimized(char*, int);
int hashFunction(int);
void insertInHashTable(bucket**, bucket*, int);
int readDemolisciStazioneParameter();
void freeBST(leaf*);
void removeFromHashTable(bucket**, int, int);
void removeFromBST(leaf**, leaf*);
leaf* nextInBST(leaf*);
leaf* minInBST(leaf*);
void readTwoIntegerParameters(int**);
void addVeichle(bucket**, int, leaf*, int);
int scrapVeichle(bucket**, int, int, int);
leaf* maxInBST(leaf*);
int getMaxAutonomyFromDistance(int, bucket**, int);
void pianificaPercorsoBackwards(leaf*, int, int**, int, bucket**, int*);
leaf* previousInBST(leaf*);
void printPianificaPercorso(int*, int);
void printIntStdoutOptimized(int);

void dijkstraForward(int, leaf*, int, bucket**, int**, int*, int*);
void enqueueWithPrio(queueElement**, queueElement*);
queueElement* removeMinFromQueue(queueElement**);
queueElement* searchInQueue(queueElement*, leaf*);
void dijkstraBackwards(int, leaf*, int, bucket**, int**, int*, int*);

// DEBUG FUNC
void inOrderBST(leaf*);
void printHashValues(bucket**);
void printQueue(queueElement*);

// const for hash values calc
const int molt = (uint32_t)((double)A * ((uint64_t)(1) << 32));

int main() {

    int end = 0;
    char* command = NULL;
    int* commandArguments = NULL;
    leaf* distanceBst = NIL;
    bucket* hashtable[HASH_DIM];
    int* output = NULL;
    int num_stations = 0;

    // initialize hash table
    for(int i = 0; i < HASH_DIM; i++) {
        hashtable[i] = NULL;
    }
    // start reading from stdin
    do{
        if(!readCommand(&command)) {    // check if the input file is empty
            #if STAMPA_STRUTTURE
                printf("BST DISTANZE:\n");
                inOrderBST(distanceBst);
                printHashValues(hashtable);
            #endif
            // TODO deallocare tutto!
            free(command);  // free all the allocated pointers
            free(commandArguments);
            end = 1;    // set flag to close software
            #if DEBUG
                printf("CHIUSURA SOFTWARE\n");
            #endif
        }
        else {
            if(stringCompare(command, AGG_STAZ, AGG_STAZ_DIM)) {    // check what type of command
                #if DEBUG
                    printf("Aggiungi stazione!\n");
                #endif
                readAggiungiStazioneParameters(&commandArguments);  // read arguments for aggiungi-stazione
                #if DEBUG
                    printf("DISTANCE: %d\n", commandArguments[0]);
                    printf("# VEICHLES: %d\n", commandArguments[1]);
                    if(commandArguments[1] > 0) {
                        for(int i = 2; i < commandArguments[1] + 2; i++) {
                            printf("VEICHLE %d AUTONOMY: %d\n", i-1, commandArguments[i]);
                        }
                    }
                #endif
                leaf* found = searchInBST(distanceBst, commandArguments[0]);    // check if already exist a station
                #if DEBUG
                    printf("END SEARCH IN BST!\n");
                    if(found == NIL) {
                        printf("FOUND = NIL\n");
                    }
                    else {
                        printf("FOUND DIST. = %d", found->key);
                    }
                #endif
                if(found == NIL) {
                    #if DEBUG
                        printf("ADDED STATION!\n");
                    #endif
                    leaf* distanceLeaf = (leaf*) malloc(sizeof(leaf));  // if not found add it to the stations BST
                    distanceLeaf->key = commandArguments[0];
                    insertInBST(&distanceBst, distanceLeaf);
                    bucket* hashTableElement = (bucket*) malloc(sizeof(bucket));    // and add his veichles into hash table
                    hashTableElement->distance = commandArguments[0];
                    hashTableElement->maxAutonomy = 0;
                    if(commandArguments[1] > 0) {
                        for(int i = 2; i < commandArguments[1] + 2; i++) {
                            leaf* veichleLeaf = (leaf*) malloc(sizeof(leaf));
                            veichleLeaf->key = commandArguments[i];
                            insertInBST(&(hashTableElement->veichles), veichleLeaf);
                            if(commandArguments[i] > hashTableElement->maxAutonomy) {
                                hashTableElement->maxAutonomy = commandArguments[i];
                            }
                        }
                    }
                    insertInHashTable(hashtable, hashTableElement, hashFunction(commandArguments[0]));
                    printStdoutOptimized(AGGIUNTA, AGGIUNTA_DIM);    // print aggiunta
                    num_stations++;
                }
                else {
                    #if DEBUG
                        printf("NOT ADDED STATION!\n");
                    #endif
                    printStdoutOptimized(NON_AGGI, NON_AGGI_DIM);    // if already exist, print non aggiunta
                }
            }
            else if(stringCompare(command, DEM_STAZ, DEM_STAZ_DIM)) {
                #if DEBUG
                    printf("Demolisci stazione!\n");
                #endif 
                int toDemolish = readDemolisciStazioneParameter();
                leaf* found = searchInBST(distanceBst, toDemolish);
                if(found != NULL) {
                    removeFromHashTable(hashtable, hashFunction(toDemolish), toDemolish);
                    removeFromBST(&distanceBst, found);
                    printStdoutOptimized(DEMOLITA, DEMOLITA_DIM);
                    num_stations--;
                }
                else {
                    printStdoutOptimized(NON_DEMO, NON_DEMO_DIM);
                }
            }
            else if(stringCompare(command, AGG_AUTO, AGG_AUTO_DIM)) {
                #if DEBUG
                    printf("Aggiungi auto!\n");
                #endif
                readTwoIntegerParameters(&commandArguments);
                leaf* found = searchInBST(distanceBst, commandArguments[0]);
                if(found != NULL) {
                    leaf* veichleToAdd = malloc(sizeof(leaf));
                    veichleToAdd->key = commandArguments[1];
                    addVeichle(hashtable, hashFunction(commandArguments[0]), veichleToAdd, commandArguments[0]);
                    printStdoutOptimized(AGGIUNTA, AGGIUNTA_DIM);
                }
                else {
                    printStdoutOptimized(NON_AGGI, NON_AGGI_DIM);
                }
            }
            else if(stringCompare(command, ROT_AUTO, ROT_AUTO_DIM)) {
                #if DEBUG
                    printf("Rottama auto!\n");
                #endif
                readTwoIntegerParameters(&commandArguments);
                leaf* found = searchInBST(distanceBst, commandArguments[0]);
                if(found != NULL) {
                    if(scrapVeichle(hashtable, hashFunction(commandArguments[0]), commandArguments[1], commandArguments[0])) {
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
                #if DEBUG
                    printf("Pianifica percorso!\n");
                #endif
                readTwoIntegerParameters(&commandArguments);
                #if DEBUG
                    printf("pianifica-percorso %d %d\n", commandArguments[0], commandArguments[1]);
                    printHashValues(hashtable);
                #endif
                if(commandArguments[0] == commandArguments[1]) {
                    printIntStdoutOptimized(commandArguments[0]);
                    putchar_unlocked('\n');
                }
                else if(commandArguments[0] < commandArguments[1]) {
                    #if DEBUG
                        printf("PIANIFICA PERCORSO ANDATA\n");
                    #endif
                    leaf* startLeaf = searchInBST(distanceBst, commandArguments[0]);
                    #if DA
                        printf("PP %d -> %d\n", commandArguments[0], commandArguments[1]);
                    #endif
                    #if DEBUG
                        printf("STAZIONE START TROVATA\n");
                    #endif
                    output = (int*) malloc(sizeof(int) * num_stations);
                    int exist = 0;
                    int count = 0;
                    #if DEBUG
                        printf("FUNZIONE PIANIFICA PERCORSO FORWARD\n");
                    #endif
                    dijkstraForward(num_stations, startLeaf, commandArguments[1], hashtable, &output, &count, &exist);
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
                    #if DEBUG
                        printf("PIANIFICA PERCORSO RITORNO\n");
                    #endif
                    leaf* startLeaf = searchInBST(distanceBst, commandArguments[0]);
                    #if DEBUG
                        printf("STAZIONE START TROVATA\n");
                    #endif
                    output = (int*) malloc(sizeof(int) * num_stations);
                    int exist = 0;
                    int count = 0;
                    #if DEBUG
                        printf("FUNZIONE PIANIFICA PERCORSO BACKWARDS\n");
                    #endif
                    dijkstraBackwards(num_stations, startLeaf, commandArguments[1], hashtable, &output, &count, &exist);
                    if(!exist) {
                        printStdoutOptimized(NESS_PER, NESS_PER_DIM);
                    }
                    else {
                        printPianificaPercorso(output, count);
                    }
                    // free the output vector mallocated
                    free(output);
                }
            }
            else {
                #if DEBUG
                    printf("Errore nella lettura di un comando!\n");
                #endif
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
                    #if MALLOC
                        printf("STO DEALLOCANDO\n");
                    #endif
                    free(*cmd);
                }
                #if MALLOC
                    printf("STO MALLOCANDO\n");
                #endif
                (*cmd) = (char*) malloc(sizeof(char) * MAX_CMD);
            }
            (*cmd)[count] = read;   // save readed char
            count++;
        }
    } while(read != ' ');   // stop when read a space
    #if DEBUG
        printf("COMANDO LETTO: %s\n", (*cmd));
    #endif
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
    #if DEBUG
        printf("DISTANCE READ!\n");
    #endif
    do{ // read from stdin # of veichles
        read = getchar_unlocked();
        if(read != ' ' && read != '\n' && read != EOF) {
            size = size*10 + read - '0';
        }
    }while(read != ' ' && read != '\n' && read != EOF);
    #if DEBUG
        printf("SIZE READ, %d!\n", size);
    #endif
    if(*vector != NULL) {
        free(*vector);
    }
    #if DEBUG
        printf("POST FREE!\n");
    #endif
    (*vector) = (int*) malloc(sizeof(int) * (size + 2));    // allocate memory
    (*vector)[0] = distance;    // save as first value # veichles
    (*vector)[1] = size;
    #if DEBUG
        printf("ALLOCATED VECTOR AND ASSIGNED VALUES!\n");
    #endif
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

    #if DEBUG
        printf("CERCANDO NEL BST LA STAZIONE CON DISTANZA: %d\n", x);
    #endif
    if(T == NIL || T->key == x) {
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

int hashFunction(int k) {   // calc of hash function

    #if DEBUG
    int h = (k * molt) % HASH_DIM;
        printf("CALCOLO HASH PER %d -> %d\n", k, h);
    #endif
    return abs(k * molt) % HASH_DIM;

}

int abs(int x) {

    return x < 0 ? (-1) * x : x;

}

void insertInHashTable(bucket** hashTable, bucket* element, int pos) {

    #if DEBUG 
        printf("AGGIUNGO ELEMENTO IN POS %d\n", pos);
    #endif
    if(hashTable[pos] != NULL) {    // if element already exist, add it to the head of the list
        #if DEBUG
            printf("ELEMENTO GIA' PRESENTE NELLA POSIZIONE, LO AGGIUNGO ALLA LISTA");
        #endif
        element->next = (struct bucket*) hashTable[pos];
    }
    hashTable[pos] = element;
    return;

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
    free(T);
    return;

}

void removeFromHashTable(bucket** hashTable, int hash, int distance) {  // remove element from hash table
    
    bucket* curr = hashTable[hash];
    if(curr->next == NULL) {    // if there is no list, free the element
        freeBST(curr->veichles);
        free(curr);
        hashTable[hash] = NULL;
    }
    else { 
        bucket* pre = NULL;
        while(curr->distance != distance) { // search for the right element in the list
            pre = curr;
            curr = (bucket*) curr->next;
        }
        if(pre != NULL) {   // if the element to free is not the head, set the list continuity
            pre->next = curr->next;
        }
        else {  // else if it's the head, set the next as head on the hash table
            hashTable[hash] = (bucket*) curr->next;
        }
        freeBST(curr->veichles);
        free(curr);
    }
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

void addVeichle(bucket** hashTable, int hash, leaf* veichleToAdd, int distance) {   // adds a new veichle into the BST of veichles of one station

    bucket* curr = hashTable[hash]; // start from the head of the list

    while(curr->distance != distance) {
        curr = (bucket*) curr->next;    // search for the bucket with the right distance
    }
    insertInBST(&(curr->veichles), veichleToAdd);   // add new veichle in BST
    if(veichleToAdd->key > curr->maxAutonomy) {
        curr->maxAutonomy = veichleToAdd->key;
    }

}

int scrapVeichle(bucket** hashTable, int hash, int autonomy, int distance) {    // remove, if exists, an element from the BST

    bucket* curr = hashTable[hash]; // statr from the head of the list

    while(curr->distance != distance) {
        curr = (bucket*) curr->next;    // search for the bucket with the right distance
    }
    leaf* toScrap = searchInBST(curr->veichles, autonomy);  // check if exists the element to scrap
    if(toScrap == NULL) {   // if not exists return false
        return 0;
    }
    else {
        int toScrapAutonomy = toScrap->key;         // save the distance
        removeFromBST(&curr->veichles, toScrap);    // else remove it and return true
        if(curr->maxAutonomy == toScrapAutonomy) {
            curr->maxAutonomy = maxInBST(curr->veichles)->key;
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

int getMaxAutonomyFromDistance(int dist, bucket** hashTable, int hash) {
    
    bucket* curr = hashTable[hash];
    while(curr->distance != dist) {
        curr = (bucket*) curr->next;
    }
    return curr->maxAutonomy;

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

void dijkstraForward(int nStations, leaf* s, int dest, bucket** hashTable, int** out, int* count, int* exist) {

    #if DA
        printf("DIJKSTRA\n");
    #endif
    queueElement* Q = NULL;
    leaf* v = s;
    while(v->key <= dest) {
        #if DA
            printf("INSERISCO IN CODA %d\n", v->key);
        #endif
        queueElement* el = (queueElement*) malloc(sizeof(queueElement));
        el->element = v;
        if(v == s) {
            el->distance = 0;
        } 
        else {
            el->distance = MAX_INT;
            el->previous = NULL;
        }
        enqueueWithPrio(&Q, el);
        if(v->key != dest) {
            v = nextInBST(v);
        }
        else {
            break;
        }
    }
    #if DA
        printf("FINE INSERIMENTO IN CODA\n");
        printQueue(Q);
    #endif
    #if DA
        printf(" - START DIJKSTRA\n");
    #endif
    while(Q != NULL) {
        queueElement* u = removeMinFromQueue(&Q);
        int uMaxAutonomy = getMaxAutonomyFromDistance(u->element->key, hashTable, hashFunction(u->element->key));
        #if DA
            printf(" * RIMOSSO MIN: %d - %d\n", u->element->key, u->distance);
            printQueue(Q);
        #endif
        v = nextInBST(u->element);
        if(v->key > u->element->key + uMaxAutonomy) {
            return;
        }
        while(v->key <= u->element->key + uMaxAutonomy && v->key <= dest) {
            #if DA
                printf("SUCCESSIVO DI %d, AUTONOMIA %d: %d\n", u->element->key, uMaxAutonomy, v->key);
            #endif
            int ndist = u->distance + 1;
            queueElement* vQEl = searchInQueue(Q, v);
            if(vQEl != NULL && vQEl->distance > ndist ) {
                vQEl->distance = ndist;
                vQEl->previous = (struct queueElement*) u;
                #if DA
                    printf("CODA AGGIORNATA!\n");
                    printQueue(Q);
                #endif
            }
            if(v->key != dest) {
                v = nextInBST(v);
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
                return;
            }
        }
    }
    #if DA
        printf("CODA VUOTA!");
    #endif

}

void enqueueWithPrio(queueElement** Q, queueElement* v)  {

    if(*Q == NULL) {
        *Q = v;
    }
    else {
        queueElement* curr = *Q;
        queueElement* prev = NULL;
        while(curr->next != NULL) {
            if(curr->distance >= v->distance) {
                prev->next = (struct queueElement*) v;
                v->next = (struct queueElement*) curr;
                return;
            }
            prev = curr;
            curr = (queueElement*) curr->next;
        }
        curr->next = (struct queueElement*) v; 
    }
    return;

}

queueElement* removeMinFromQueue(queueElement** Q) {

    queueElement* curr = *Q;
    queueElement* pred = NULL;
    queueElement* min = NULL;
    queueElement* predMin = NULL;

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
        else if(curr->distance == min->distance && curr->element->key < min->element->key) {
            min = curr;
            predMin = pred;
        }
        pred = curr;
        curr = (queueElement*) curr->next;
    }
    if(predMin == NULL) {
        *Q = (queueElement*) min->next;
    }
    else {
        predMin->next = min->next;
    }
    return min;

}

queueElement* searchInQueue(queueElement* Q, leaf* q) {

    queueElement* tmp = Q;
    while(tmp->element != q) {
        if(tmp->element != q && tmp->next == NULL) {
            return NULL;
        }
        tmp = (queueElement*) tmp->next;
    }
    return tmp;

}

void dijkstraBackwards(int nStations, leaf* s, int dest, bucket** hashTable, int** out, int* count, int* exist) {

    #if DA
        printf("DIJKSTRA\n");
    #endif
    queueElement* Q = NULL;
    leaf* v = s;
    while(v->key >= dest) {
        #if DA
            printf("INSERISCO IN CODA %d\n", v->key);
        #endif
        queueElement* el = (queueElement*) malloc(sizeof(queueElement));
        el->element = v;
        if(v == s) {
            el->distance = 0;
        } 
        else {
            el->distance = MAX_INT;
            el->previous = NULL;
        }
        enqueueWithPrio(&Q, el);
        if(v->key != dest) {
            v = previousInBST(v);
        }
        else {
            break;
        }
    }
    #if DA
        printf("FINE INSERIMENTO IN CODA\n");
        printQueue(Q);
    #endif
    #if DA
        printf(" - START DIJKSTRA\n");
    #endif
    while(Q != NULL) {
        queueElement* u = removeMinFromQueue(&Q);
        int uMaxAutonomy = getMaxAutonomyFromDistance(u->element->key, hashTable, hashFunction(u->element->key));
        #if DA
            printf(" * RIMOSSO MIN: %d - %d\n", u->element->key, u->distance);
            printQueue(Q);
        #endif
        v = previousInBST(u->element);
        if(v->key < u->element->key - uMaxAutonomy) {
            return;
        }
        while(v->key >= u->element->key - uMaxAutonomy && v->key >= dest) {
            #if DA
                printf("SUCCESSIVO DI %d, AUTONOMIA %d: %d\n", u->element->key, uMaxAutonomy, v->key);
            #endif
            int ndist = u->distance + 1;
            queueElement* vQEl = searchInQueue(Q, v);
            if(vQEl != NULL && vQEl->distance > ndist ) {
                vQEl->distance = ndist;
                vQEl->previous = (struct queueElement*) u;
                #if DA
                    printf("CODA AGGIORNATA!\n");
                    printQueue(Q);
                #endif
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
                return;
            }
        }
    }
    #if DA
        printf("CODA VUOTA!");
    #endif  

}

// ********* DEBUG FUNCTIONS *********

void inOrderBST(leaf* T) {

    if(T->left != NIL) {
        inOrderBST((leaf*) T->left);
    }
    printf("%d\n", T->key);
    if(T->right != NIL) {
        inOrderBST((leaf*) T->right);
    }
    return;

}

void printHashValues(bucket** hashTable) {

    for(int i = 0; i < HASH_DIM; i++) {
        if(hashTable[i] != NULL) {
            bucket* curr = hashTable[i];
            while(curr != NULL) {
                // printf("ELEMENTO IN POS %d\n", i);
                printf("%d", curr->distance);
                // printf("PARCO AUTO:\n");
                // inOrderBST(curr->veichles);
                printf(" - MAX_A: %d\n", curr->maxAutonomy);
                curr = (bucket*) curr->next;
            }
        }
    }
    return;

}

void printQueue(queueElement* Q) {

    queueElement* tmp = Q;
    while(tmp != NULL) {
        printf("ELEMENTO %d - DISTANCE %d\n", tmp->element->key, tmp->distance);
        tmp = (queueElement*) tmp->next;
    }
    return;

}