#include "GridTools.h"
double inf = numeric_limits<double>::infinity();
double heuristicVal(imat node, int h){

        // number of empty cells
        if(h == 1){
            imat indexes = indexOfElement(node, 0);
            int numOfBlank = indexes.n_rows;
            return numOfBlank;
        }

        // monotonicity
        else if (h == 2) {
            int rowNum = 0, colNum = 0;
            double orderScore = 0;
            for (rowNum = 0; rowNum < 4; rowNum++) {
                orderScore += getMonotonicityOfVector(node.row(rowNum).t());
                if (rowNum == 0 || rowNum == 3) {
                    orderScore += getMonotonicityOfVector(node.row(rowNum).t());
                }
            }
            for (colNum = 0; colNum < 4; colNum++) {
                orderScore += getMonotonicityOfVector(node.col(colNum));
                if (colNum == 0 || colNum == 3) {
                    orderScore += getMonotonicityOfVector(node.col(colNum));
                }
            }
            return orderScore;
        }

        // constant Snake with all other heuristics
        else if(h == 3){
            imat indexes = indexOfElement(node, 0);
            int numOfBlank = indexes.n_rows;
            mat weight1 = {{pow(10,4), pow(10,3), pow(10,2), pow(10,1)},
                           {pow(10,5), pow(10,6), pow(10,7), pow(10,8)},
                           {pow(10,12), pow(10,11), pow(10,10), pow(10,9)},
                           {pow(10,13), pow(10,14), pow(10,15), pow(10,16)}};

            mat weight2 = {{pow(10,4), pow(10,5), pow(10,12), pow(10,13)},
                           {pow(10,3), pow(10,6), pow(10,11), pow(10,14)},
                           {pow(10,2), pow(10,7), pow(10,10), pow(10,15)},
                           {pow(10,1), pow(10,8), pow(10,9), pow(10,16)}};

            mat weight3 = {{pow(10,1), pow(10,2), pow(10,3), pow(10,4)},
                           {pow(10,8), pow(10,7), pow(10,6), pow(10,5)},
                           {pow(10,9), pow(10,10), pow(10,11), pow(10,12)},
                           {pow(10,16), pow(10,15), pow(10,14), pow(10,13)}};

            mat weight4 = {{pow(10,1), pow(10,8), pow(10,9), pow(10,16)},
                           {pow(10,2), pow(10,7), pow(10,10), pow(10,15)},
                           {pow(10,3), pow(10,6), pow(10,11), pow(10,14)},
                           {pow(10,4), pow(10,5), pow(10,12), pow(10,13)}};

            double snake [] = {0, 0, 0, 0};
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    snake[0] += node(i, j) * weight1(i, j);
                    snake[1] += node(i, j) * weight2(i, j);
                    snake[2] += node(i, j) * weight3(i, j);
                    snake[3] += node(i, j) * weight4(i, j);
                }
            }
            double sumOfPairwise = 0;
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++) {
                    sumOfPairwise += abs(node(i,j)-node(i,j+1));
                    sumOfPairwise += abs(node(i,j)-node(i+1,j));
                }
            }
            for(int i = 1; i < 4; i++){
                if(snake[i] > snake[0])
                    snake[0] = snake[i];
            }
            ivec temp(vectorise(node));
            temp = sort(temp);
            double penalty [] = {0, 0};
            if (node[3,3] != temp[15])
                penalty[0] = temp[15] + temp[14] + temp[13] + temp[12] - node[3, 3];
            else if(node[3,2] != temp[14])
                penalty[0] = temp[14] + temp[13] + temp[12] - node[3, 2];
            else if(node[3,1] != temp[13])
                penalty[0] = temp[13] + temp[12] - node[3, 1];
            else if(node[3,0] != temp[12])
                penalty[0] = temp[12] - node[3, 0];

            if (node[3, 3] != temp[15])
                penalty[1] = temp[15] + temp[14] + temp[13] + temp[12] - node[3, 3];
            else if(node[2, 3] != temp[14])
                penalty[1] = temp[14] + temp[13] + temp[12] - node[2, 3];
            else if(node[1, 3] != temp[13])
                penalty[1] = temp[13] + temp[12] - node[1, 3];
            else if (node[0, 3] != temp[12])
                penalty[1] = temp[12] - node[0, 3];

            if(penalty[1] < penalty[0])
                penalty[0] = penalty[1];

            if (node[3, 0] != temp[15])
                penalty[1] = temp[15] + temp[14] + temp[13] + temp[12] - node[3, 0];
            else if(node[3, 1] != temp[14])
                penalty[1] = temp[14] + temp[13] + temp[12] - node[3, 1];
            else if(node[3, 2] != temp[13])
                penalty[1] = temp[13] + temp[12] - node[3, 2];
            else if (node[3, 3] != temp[12])
                penalty[1] = temp[12] - node[3, 3];

            if(penalty[1] < penalty[0])
                penalty[0] = penalty[1];

            if (node[0, 3] != temp[15])
                penalty[1] = temp[15] + temp[14] + temp[13] + temp[12] - node[0, 3];
            else if(node[1, 3] != temp[14])
                penalty[1] = temp[14] + temp[13] + temp[12] - node[1, 3];
            else if(node[2, 3] != temp[13])
                penalty[1] = temp[13] + temp[12] - node[2, 3];
            else if (node[3, 3] != temp[12])
                penalty[1] = temp[12] - node[3, 3];

            if(penalty[1] < penalty[0])
                penalty[0] = penalty[1];

            if(numOfBlank == 0)
                return -inf;
            return pow(snake[0],3) - pow(sumOfPairwise, 3)- pow(penalty[0],3);

        }

        // dynamic snakes
        else if (h == 4) {

            int base = 2, snakeNum = 0, snakeValue = 0, bestSnakeValue = 0;
            imat snakePit [8] = {   imat { {pow(base, 4), pow(base, 3), pow(base, 2), 1},
                                           {pow(base, 5), pow(base, 6), pow(base, 7), pow(base, 8)},
                                           {pow(base, 12), pow(base, 11), pow(base, 10), pow(base, 9)},
                                           {pow(base, 13), pow(base, 14), pow(base, 15), pow(base, 16)} },

                                    imat { {1, pow(base, 2), pow(base, 3), pow(base, 4)},
                                           {pow(base, 8), pow(base, 7), pow(base, 6), pow(base, 5)},
                                           {pow(base, 9), pow(base, 10), pow(base, 11), pow(base, 12)},
                                           {pow(base, 16), pow(base, 15), pow(base, 14), pow(base, 13)} },

                                    imat { {pow(base, 16), pow(base, 15), pow(base, 14), pow(base, 13)},
                                           {pow(base, 9), pow(base, 10), pow(base, 11), pow(base, 12)},
                                           {pow(base, 8), pow(base, 7), pow(base, 6), pow(base, 5)},
                                           {1, pow(base, 2), pow(base, 3), pow(base, 4)} },

                                    imat { {pow(base, 13), pow(base, 14), pow(base, 15), pow(base, 16)},
                                           {pow(base, 12), pow(base, 11), pow(base, 10), pow(base, 9)},
                                           {pow(base, 5), pow(base, 6), pow(base, 7), pow(base, 8)},
                                           {pow(base, 4), pow(base, 3), pow(base, 2), 1} },

                                    imat { {pow(base, 4), pow(base, 5), pow(base, 12), pow(base, 13)},
                                           {pow(base, 3), pow(base, 6), pow(base, 11), pow(base, 14)},
                                           {pow(base, 2), pow(base, 7), pow(base, 10), pow(base, 15)},
                                           {1, pow(base, 8), pow(base, 9), pow(base, 16)} },

                                    imat { {1, pow(base, 8), pow(base, 9), pow(base, 16)},
                                           {pow(base, 2), pow(base, 7), pow(base, 10), pow(base, 15)},
                                           {pow(base, 3), pow(base, 6), pow(base, 11), pow(base, 14)},
                                           {pow(base, 4), pow(base, 5), pow(base, 12), pow(base, 13)} },

                                    imat { {pow(base, 16), pow(base, 9), pow(base, 8), 1},
                                           {pow(base, 15), pow(base, 10), pow(base, 7), pow(base, 2)},
                                           {pow(base, 14), pow(base, 11), pow(base, 6), pow(base, 3)},
                                           {pow(base, 13), pow(base, 12), pow(base, 5), pow(base, 4)} },

                                    imat { {pow(base, 13), pow(base, 12), pow(base, 5), pow(base, 4)},
                                           {pow(base, 14), pow(base, 11), pow(base, 6), pow(base, 3)},
                                           {pow(base, 15), pow(base, 10), pow(base, 7), pow(base, 2)},
                                           {pow(base, 16), pow(base, 9), pow(base, 8), 1} }
            };

            for (snakeNum = 0; snakeNum < 8; snakeNum++) {
                snakeValue = sum(sum(node % snakePit[snakeNum]));
                if (snakeValue > bestSnakeValue) {
                    bestSnakeValue = snakeValue;
                }
            }

            int rowNum = 0, colNum = 0;
            double orderScore = 0;

            for (rowNum = 0; rowNum < 4; rowNum++) {
                orderScore += getMonotonicityOfVector(node.row(rowNum).t());
            }
            for (colNum = 0; colNum < 4; colNum++) {
                orderScore += getMonotonicityOfVector(node.col(colNum));
            }

            imat indexes = indexOfElement(node, 0);
            int numOfBlank = indexes.n_rows;

            return bestSnakeValue + orderScore + 30 * numOfBlank;

        }
        else {
            return 0;
        }
}

double expectimax(imat node, int ply, int h, string expOrMax){
    if(!checkIfCanGo(node))
        return -inf;
    else if(ply == 0)
        return heuristicVal(node, h);
    else if(!expOrMax.compare("max")){
      string moves[] = {"Up", "Left", "Down", "Right"};
      double maxVal = -inf;
      for(int i = 0; i < 4; i++){
        imat child = getChild(node, moves[i], false);
        double val = 0;
        if(isEqual(child, node))
          continue;
        else
          val = expectimax(child, ply - 1, h, "exp");

        if(val> maxVal)
          maxVal = val;
      }
      return maxVal;
    }

    else if(!expOrMax.compare("exp")){
        imat indexes = indexOfElement(node, 0);
        int numOfRchild = indexes.n_rows;
        double sumOf2 = 0;
        for(int i = 0; i < numOfRchild; i++){
            imat rChild = node;
            rChild(indexes(i,0), indexes(i,1)) = 2;
            sumOf2 += expectimax(rChild, ply, h, "max");
        }
        double sumOf4 = 0;
        for(int i = 0; i < numOfRchild; i++){
            imat rChild = node;
            rChild(indexes(i,0), indexes(i,1)) = 4;
            sumOf4 += expectimax(rChild, ply, h, "max");
        }
        return 0.9*(sumOf2/numOfRchild) + 0.1*(sumOf4/numOfRchild);
    }

  }
string getBestAction(imat node, int ply, int h){
    string bestMove = "";
    string moves[] = {"Up", "Left", "Down", "Right"};
    double maxVal = -inf;
    for(int i = 0; i < 4; i++){
        imat child = getChild(node, moves[i], false);
        double val = 0;
        if(isEqual(child, node))
        continue;
        if(!checkIfCanGo(child))
            val = -inf;
        else{
            val = expectimax(child, ply - 1, h, "exp");
        }
        if(val>= maxVal){
            maxVal = val;
            bestMove = moves[i];
        }

      }

      return bestMove;

  }
