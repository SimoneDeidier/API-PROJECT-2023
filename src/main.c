#include <stdlib.h>
#include <stdio.h>
#define MAX_COMMAND_DIM 19

int readCommand(char*);
int stringCompare(char*, char*, int);

int main() {

    char command[MAX_COMMAND_DIM] = {'\0'};
    int closeApp = 0;

    do {
        closeApp = readCommand(command);
        printf("%s", command);
    } while(!closeApp);

    return 0;
}

int readCommand(char* command) {
    char read = '\0';
    int count = 0;

    do {
        read = getchar_unlocked();
        if(read == '\0') {
            return 1;
        }
        if(read != ' ') {
            command[count] = read;
            count++;
        }
    } while(read != ' ');
    return 1;
}

int stringCompare(char* a, char* b, int dim) {
    for(int i = 0; i < dim; i++) {
        if(a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}