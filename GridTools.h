#include <iostream>
#include <armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string>
#include <math.h>       /* pow */
#include <limits>
using namespace std;
using namespace arma;
imat indexOfElement(imat m, int k){
    uvec q2 = find(m == k);
    imat indexes(q2.n_elem, 2, fill::zeros);
    for(int i = 0; i < q2.n_elem; i++){
      indexes(i,1) = floor(q2(i)/4);
      indexes(i,0) = q2(i)%4;
    }

    return indexes;
}
imat addRandomTiles(imat tiles){
    if (tiles.min() != 0)
      return tiles;
    imat indexes = indexOfElement(tiles, 0);
    int randomIndex = rand() % indexes.n_rows;
    int randTiles = 2;
    double r = ((double) rand() / (RAND_MAX));
    if (r <= 0.1)
      randTiles = 4;
    tiles(indexes(randomIndex,0), indexes(randomIndex,1)) = randTiles;
    return tiles;
}
imat moveTiles(imat tiles,  string move){
  if (!move.compare("Up")){
    for(int j = 0; j < 4; j++){
      for(int i = 1; i < 4; i++){
        if(tiles(i, j) == 0)
          continue;
        int temp_i = i;
        while(temp_i > 0 ){
          if(tiles(temp_i - 1, j) != 0)
            break;
          int temp = tiles(temp_i, j);
          tiles(temp_i, j) = tiles(temp_i - 1, j);
          tiles(temp_i - 1, j) = temp;
          temp_i -= 1;
        }
      }
    }
  }
  else if(!move.compare("Left")){
    for(int i = 0; i < 4; i++){
      for(int j = 1; j < 4; j++){
        if(tiles(i, j) == 0)
          continue;
        int temp_j = j;
        while(temp_j > 0 ){
          if(tiles(i, temp_j - 1) != 0)
            break;
          int temp = tiles(i, temp_j);
          tiles(i, temp_j) = tiles(i, temp_j - 1);
          tiles(i, temp_j - 1) = temp;
          temp_j -= 1;
        }

      }
    }

  }
  else if(!move.compare("Down")){
    for(int j = 0; j < 4; j++){
      for(int i = 2; i >= 0; i--){
        if(tiles(i, j) == 0)
          continue;
        int temp_i = i;
        while(temp_i <= 2 ){
          if(tiles(temp_i + 1, j) != 0)
            break;
          int temp = tiles(temp_i, j);
          tiles(temp_i, j) = tiles(temp_i + 1, j);
          tiles(temp_i + 1, j) = temp;
          temp_i += 1;
        }
      }
    }
  }
  else if(!move.compare("Right")){
    for(int i = 0; i < 4; i++){
      for(int j = 2; j >= 0; j--){
        if(tiles(i, j) == 0)
          continue;
        int temp_j = j;
        while(temp_j <= 2 ){
          if(tiles(i, temp_j + 1) != 0)
            break;
          int temp = tiles(i, temp_j);
          tiles(i, temp_j) = tiles(i, temp_j + 1);
          tiles(i, temp_j + 1) = temp;
          temp_j += 1;
        }

      }
    }

  }
  return tiles;
}
imat mergeTiles(imat tiles,string move){
  if (!move.compare("Up")){
    for(int j = 0; j < 4; j++){
      for(int i = 0; i < 3; i++){
        if(tiles(i, j) == 0 || tiles(i, j) != tiles(i+1, j))
          continue;
        tiles(i, j) *= 2;
        tiles(i+1, j) = 0;
      }
    }
  }
  else if(!move.compare("Left")){
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 3; j++){
        if(tiles(i, j) == 0 || tiles(i, j) != tiles(i, j+1))
          continue;
        tiles(i, j) *= 2;
        tiles(i, j+1) = 0;
      }
    }

  }
  else if(!move.compare("Down")){
    for(int j = 0; j < 4; j++){
      for(int i = 3; i > 0; i--){
        if(tiles(i, j) == 0 || tiles(i, j) != tiles(i-1, j))
          continue;
        tiles(i, j) *= 2;
        tiles(i-1, j) = 0;
      }
    }
  }
  else if(!move.compare("Right")){
    for(int i = 0; i < 4; i++){
      for(int j = 3; j > 0; j--){
        if(tiles(i, j) == 0 || tiles(i, j) != tiles(i, j-1))
          continue;
        tiles(i, j) *= 2;
        tiles(i, j-1) = 0;
      }
    }

  }
  tiles = moveTiles(tiles, move);
  return tiles;
}
bool isEqual(imat a, imat b){

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      if(a(i, j) != b(i, j))
        return false;
    }
  }
  return true;
}
bool canMove(imat tiles, string move){
  imat temp = tiles;
  tiles = moveTiles(tiles, move);
  tiles = mergeTiles(tiles, move);

  if(! isEqual(temp, tiles)){
    return true;
  }

  return false;
}
bool checkIfCanGo(imat tiles){
  if(tiles.min() != 0){
    string move[] = {"Up","Right", "Left", "Down"};
    for(int i = 0; i < 4; i++){
      if (canMove(tiles, move[i]))
        return true;
    }
    return false;
  }
  return true;
}
imat getChild(imat tiles, string move, bool addrandom){

    if(canMove(tiles, move)){
      tiles = moveTiles(tiles, move);
      tiles = mergeTiles(tiles, move);
      if(addrandom)
        tiles = addRandomTiles(tiles);
    }
    return tiles;
}
double getMonotonicityOfVector(ivec vector) {
  bool dec = true, inc = true;
  int elNum = 0;
  double dist = 0;

  for (elNum = 0; elNum < 3; elNum++) {
    if (vector(elNum) == 0) {
      vector(elNum) = 1;
    }
    if (vector(elNum + 1) == 0) {
      vector(elNum + 1) = 1;
    }
    if (dec && vector(elNum) < vector(elNum + 1)) {
      dec = false;
    }
    else if (inc && vector(elNum) > vector(elNum + 1)) {
      inc = false;
    }
    dist += abs(log2(vector(elNum)) - log2(vector(elNum + 1)));
  }

  if (dec && inc) {
    return 0;
  }
  else if (dec || inc) {
    return sum(vector) / dist;
  }
  else {
    return -1;
  }
}
