#include "AI.h"
void printTiles(imat tiles){
        cout<<"****************************************************************************\n";
        cout<<"****************************************************************************\n";
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++) {
                if (tiles(i, j) == 0) {
                    cout << "_\t\t\t";
                }
                else {
                    cout << tiles(i, j) << "\t\t\t";
                }
            }
            cout<<"\n";
        }
}
int main(){
        int ply, h;
        cout << "Enter number of ply: ";
        cin >> ply;
        cout<<"Enter heuristic number(1 to 4): ";
        cin >>h;
        srand(time(NULL));
        int i = 0;
        imat tiles(4, 4, fill::zeros);
        tiles = addRandomTiles(tiles);
        tiles = addRandomTiles(tiles);
        printTiles(tiles);
        string move = "";
        while(checkIfCanGo(tiles)){
            move = getBestAction(tiles, ply, h);
            tiles = getChild(tiles, move, true);
            printTiles(tiles);
        }
        cout<<"GAME OVER\n";
        cout<<"Max tile value is "<<max(max(tiles))<<".";




}
