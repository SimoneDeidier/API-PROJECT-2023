#include <stdlib.h>
#include <stdio.h>

#define MAX_CMD 19
#define AGG_STAZ "aggiungi-stazione"
#define DEM_STAZ "demolisci-stazione"
#define AGG_AUTO "aggiungi-auto"
#define ROT_AUTO "rottama-auto"
#define PIA_PERC "pianifica-percorso"
#define AGG_STAZ_DIM 17
#define DEM_STAZ_DIM 18
#define AGG_AUTO_DIM 13
#define ROT_AUTO_DIM 12
#define PIA_PERC_DIM 18
#define NIL NULL

#define DEBUG 1
#define MALLOC 0

typedef struct {

    struct leaf* left;
    struct leaf* right;
    struct leaf* p;
    int key;

}leaf;

int readCommand(char**);
int stringCompare(char*, char*, int);
void readAggiungiStazioneParameters(int**);
void insertInBST(leaf**, leaf**);

int main() {

    int end = 0;
    //int* arguments = NULL;
    char* command = NULL;
    int* commandArguments = NULL;
    leaf* distanceBst = NIL;
    do{
        if(!readCommand(&command)) {    // check if the input file is empty
            #if DEBUG
                printf("CHIUSURA SOFTWARE\n");
            #endif
            // TODO deallocare tutto!
            free(command);  // free all the allocated pointers
            free(commandArguments);
            end = 1;    // set flag to close software
        }
        else {
            if(stringCompare(command, AGG_STAZ, AGG_STAZ_DIM)) {    // check what type of command
                #if DEBUG
                    printf("Aggiungi stazione!\n");
                #endif
                readAggiungiStazioneParameters(&commandArguments);
                #if DEBUG
                    printf("DISTANCE: %d\n", commandArguments[0]);
                    printf("# VEICHLES: %d\n", commandArguments[1]);
                    for(int i = 2; i < commandArguments[1] + 2; i++) {
                        printf("VEICHLE %d AUTONOMY: %d\n", i, commandArguments[i]);
                    }
                #endif
                leaf* newLeaf = (leaf*) malloc(sizeof(leaf));
                newLeaf->key = commandArguments[0];
                insertInBST(&distanceBst, &newLeaf);
            }
            else if(stringCompare(command, DEM_STAZ, DEM_STAZ_DIM)) {
                #if DEBUG
                    printf("Demolisci stazione!\n");
                #endif  
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
                #if DEBUG && MALLOC
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
    // TODO qua si può migliorare qualcosa???
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
    (*vector) = (int*) malloc(sizeof(int) * (size + 2));    // allocate memory
    (*vector)[0] = distance;    // save as first value # veichles
    (*vector)[1] = size;
    for(int i = 2; i < size + 2; i++) { // read all veichle values from stdin
        read = '\0';
        val = 0;
        do {
            read = getchar_unlocked();
            if(read != ' ') {
                val = val*10 + read - '0';
            }
        }while(read != ' ');
        (*vector)[i] = val;
    }
}

void insertInBST(leaf** T, leaf** x) {
    leaf* pre = NIL;
    leaf* cur = (*T);

    while(cur != NIL) {
        pre = cur;
        if((*x)->key < cur->key) {
            cur = (leaf*) cur->left;
        }
        else {
            cur = (leaf*) cur->right;
        }
    }
    (*x)->p = (struct leaf*) pre;
    if(pre == NIL) {
        (*T) = (*x);
    }
    else if((*x)->key < pre->key) {
        pre->left = (struct leaf*) (*x);
    }
    else {
        pre->right = (struct leaf*) (*x);
    }
}