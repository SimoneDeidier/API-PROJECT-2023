#include <stdlib.h>
#include <stdio.h>

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
#define AGGIUNTA_DIM 8
#define NON_AGGI_DIM 12
#define DEMOLITA_DIM 8
#define NON_DEMO_DIM 12

#define HASH_DIM 1024 * 16
#define A 0.6180339887498948482045868343656381177203091798057628621354486227

#define DEBUG 0
#define MALLOC 0
#define STAMPA_STRUTTURE 1

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

}bucket;

int readCommand(char**);
int stringCompare(char*, char*, int);
void readAggiungiStazioneParameters(int**);
void insertInBST(leaf**, leaf*);
leaf* searchInBST(leaf*, int);
void printStdinOptimized(char*, int);
int hashFunction(int);
void insertInHashTable(bucket**, bucket*, int);
int readDemolisciStazioneParameter();
void freeBST(leaf*);
void removeFromHashTable(bucket**, int, int);

// DEBUG FUNC
void inOrderBST(leaf*);
void printHashValues(bucket**);

// const for hash values calc
const int molt = (uint32_t)((double)A * ((uint64_t)1 << 32));

int main() {

    int end = 0;
    char* command = NULL;
    int* commandArguments = NULL;
    leaf* distanceBst = NIL;
    bucket* hashtable[HASH_DIM];

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
                    for(int i = 2; i < commandArguments[1] + 2; i++) {
                        printf("VEICHLE %d AUTONOMY: %d\n", i-1, commandArguments[i]);
                    }
                #endif
                leaf* found = searchInBST(distanceBst, commandArguments[0]);    // check if already exist a station
                #if DEBUG
                    printf("END SEARCH IN BST!\n");
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
                    for(int i = 2; i < commandArguments[1] + 2; i++) {
                        leaf* veichleLeaf = (leaf*) malloc(sizeof(leaf));
                        veichleLeaf->key = commandArguments[i];
                        insertInBST(&(hashTableElement->veichles), veichleLeaf);
                    }
                    insertInHashTable(hashtable, hashTableElement, hashFunction(commandArguments[0]));
                    printStdinOptimized(AGGIUNTA, AGGIUNTA_DIM);    // print aggiunta
                }
                else {
                    #if DEBUG
                        printf("NOT ADDED STATION!\n");
                    #endif
                    printStdinOptimized(NON_AGGI, NON_AGGI_DIM);    // if already exist, print non aggiunta
                }
            }
            else if(stringCompare(command, DEM_STAZ, DEM_STAZ_DIM)) {
                #if DEBUG
                    printf("Demolisci stazione!\n");
                #endif 
                int toDemolish = readDemolisciStazioneParameter();
                leaf* found = searchInBST(distanceBst, toDemolish);
                if(found != NULL) {
                    // TODO rimuovere da bst la distanza
                    removeFromHashTable(hashtable, hashFunction(toDemolish), toDemolish);
                    printStdinOptimized(DEMOLITA, DEMOLITA_DIM);
                }
                else {
                    printStdinOptimized(NON_DEMO, NON_DEMO_DIM);
                }
            }
            else if(stringCompare(command, AGG_AUTO, AGG_AUTO_DIM)) {
                #if DEBUG
                    printf("Aggiungi auto!\n");
                #endif
            }
            else if(stringCompare(command, ROT_AUTO, ROT_AUTO_DIM)) {
                #if DEBUG
                    printf("Rottama auto!\n");
                #endif
            }
            // TODO questo else if si può togliere secondo me
            else if(stringCompare(command, PIA_PERC, PIA_PERC_DIM)) {
                #if DEBUG
                    printf("Pianifica percorso!\n");
                #endif
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

    
    do{ // read from stdin # of veichles
        read = getchar_unlocked();
        if(read != ' ') {
            distance = distance*10 + read - '0';
        }
    }while(read != ' ');
    do{ // read from stdin # of veichles
        read = getchar_unlocked();
        if(read != ' ') {
            size = size*10 + read - '0';
        }
    }while(read != ' ');
    if(*vector != NULL) {
        free(*vector);
    }
    (*vector) = (int*) malloc(sizeof(int) * (size + 2));    // allocate memory
    (*vector)[0] = distance;    // save as first value # veichles
    (*vector)[1] = size;
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

void printStdinOptimized(char* str, int dim) {
    
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

void removeFromHashTable(bucket** hashTable, int hash, int distance) {
    
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
            printf("ELEMENTO IN POS %d\n", i);
            printf("PARCO AUTO:\n");
            inOrderBST(hashTable[i]->veichles);
        }
    }
    return;

}