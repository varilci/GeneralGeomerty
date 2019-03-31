#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/adjacency_list.h>
#include <igl/exact_geodesic.h>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include<bits/stdc++.h>
#include <algorithm>
#include <igl/jet.h>
#include <igl/edges.h>
#include <igl/vertex_components.h>

#include <igl/isolines.h>

#include <igl/project_to_line.h>

#include <math.h>

#include <fstream>
#include <string>


using namespace std;
# define INF 0x3f3f3f3f

Eigen::MatrixXd V;
Eigen::MatrixXi F;

typedef pair<float, int> nodeDistPair;
typedef pair<double, int> nodeDistPairDouble;


float distance(float x1, float y1, float z1, float x2, float y2, float z2) {

    return std::sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));

}

// Prints shortest path costs from src to all other vertices
vector<float> shortestPath(int src, int noOfVertices, std::vector<std::vector<double>> adjList)
{
    priority_queue< nodeDistPair, vector <nodeDistPair> , greater<nodeDistPair> > pq;

    std::vector<float> dist(noOfVertices, INF); //TODO: # OF VERTICES HERE

    pq.push(make_pair(0, src));
    dist[src] = 0;

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();

        list< pair<int, int> >::iterator i;

        for (auto i: adjList[u]) {

            int v = i;
            float cost = distance(V(u,0), V(u,1), V(u,2), V(v,0), V(v,1), V(v,2));

            if(dist[v] > dist[u] + cost) {

                dist[v] = dist[u] + cost;
                pq.push(make_pair(dist[v], v));

            }

        }
    }
    return dist;
}

vector<int> shortestPathPrint(int src, int dest, int noOfVertices, std::vector<std::vector<double>> adjList)
{
    priority_queue< nodeDistPair, vector <nodeDistPair> , greater<nodeDistPair> > pq;

    std::vector<float> dist(noOfVertices, INF); //TODO: # OF VERTICES HERE

    pq.push(make_pair(0, src));
    dist[src] = 0;

    std::vector<int> parent(noOfVertices, -1);

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();

        list< pair<int, int> >::iterator i;

        for (auto i: adjList[u]) {

            int v = i;
            float cost = distance(V(u,0), V(u,1), V(u,2), V(v,0), V(v,1), V(v,2));

            if(dist[v] > dist[u] + cost) {

                //cout << " Dijkstra u parent : " << u << "out of v " << v  << endl;

                dist[v] = dist[u] + cost;
                pq.push(make_pair(dist[v], v));

                parent[v] = u;

                if(v == dest) {
                    return parent;
                }



            }

        }
    }
    return parent;
}

//pushes path values to &path
void printPath(vector<int> parent, int j, std::vector<int>& path)
{

    // Base Case : If j is source
    if (parent[j] == - 1)
        return;

    printPath(parent, parent[j], path);


    path.push_back(j);
    //printf("%d ", j);
}

int main(int argc, char *argv[])
{
    std::string file_name;
    if (argc > 1) {

        file_name = argv[1];

        if(file_name.back() == 'f') {
            igl::readOFF( file_name, V, F);
        }
        else if (file_name.back() == 'j') {
            igl::readOBJ( file_name, V, F);
        }
    }
    else
    {
        std::cout << "Please enter file name as argument before running. (Also only .off and .obj are supported)" << std::endl;
        return 0;
    }
    


    //igl::readOFF("../man0.off", V, F);

    int program;
    while(1){

        std::cout << "Welcome to Batuhan's Program! Please select what you want to do by typing the number;" << std::endl;
        std::cout << "(1) Geodesic distance matrix for whole mesh" << std::endl;
        std::cout << "(2) Geodesic distance for two query points (Needs additional input)" << std::endl;
        std::cout << "(3) K-Iso curves & Coloring" << std::endl;
        std::cout << "(4) Bilateral Mapping ( PARTIAL)" << std::endl;
        std::cout << "(5) EXIT" << std::endl;

        cin >> program;
        if (program == 1) {
            cout << "Applying Dijkstra" << endl;

            std::vector<std::vector<double>> adj_list;
            //ROW i is the number of the vertice, and the elements of that row indicates the neighbors.
            igl::adjacency_list(F,adj_list);

            std::vector<std::vector<float>> dijkstraMatrix;

            for(int i = 0; i < V.rows(); i++)
            {
                //cout << " Dijkstra number : " << i << "out of " << V.rows()  << endl;

                    //APPLY DIJKSTRA

                    std::vector<float> result = shortestPath(i, V.rows(), adj_list);

                    dijkstraMatrix.push_back(result);

            }
            cout << "Writing matrix to output file 'outputMatrix.txt'" << endl;
            std::ofstream output("outputMatrix.txt");
            for (int i=0;i<V.rows();i++)
            {
            	for (int j=0;j<V.rows();j++)
            	{
            		output << dijkstraMatrix[i][j] << " "; // behaves like cout - cout is also a stream
            	}
            	output << "\n";
            }
            cout << "DONE" << endl;
            cout << "------------------------------" << endl;

        }
        else if (program == 2) {
            int startIndex;
            int endIndex;

            std::cout << "Please enter the index # of the starting vertice out of " << V.rows() << " vertex" << std::endl;
            cin >> startIndex;
            std::cout << "Please enter the index # of the ending vertice out of " << V.rows() << " vertex" << std::endl;
            cin >> endIndex;
            //TODO:HIGHLIGHT PATH, ALSO PRINT PATH

            Eigen::MatrixXd C; //COLOR MATRIX

            C = Eigen::MatrixXd::Constant(V.rows(),3,1); //WHITE ALL

            cout << "Applying Dijkstra" << endl;
            std::vector<std::vector<double>> adj_list;
            //ROW i is the number of the vertice, and the elements of that row indicates the neighbors.
            igl::adjacency_list(F,adj_list);

            std::vector<int >path = shortestPathPrint(startIndex, endIndex, V.rows(), adj_list);

        //    std::cout << "####" << std::endl;
        //    for (int i=0;i < V.rows();i++)
        //        std::cout << path[i] << ' ';
        //    std::cout << "####" << std::endl;


            cout << "Coloring" << endl;
            std::vector<int> pathPrint;
            pathPrint.push_back(startIndex);
            printPath(path, endIndex, pathPrint);


            Eigen::MatrixXd isoV = Eigen::MatrixXd::Zero(pathPrint.size(),3);
            Eigen::MatrixXi isoE = Eigen::MatrixXi::Zero(pathPrint.size()-1,2);

            //cout << "size " << pathPrint.size()+1 << endl;

            // std::cout << "####" << std::endl;
            // for (auto iter: pathPrint)
            //     std::cout << iter << ' ';
            // std::cout << "####" << std::endl;

            // isoV(0,0) = V(startIndex,0);
            // isoV(0,1) = V(startIndex,1);
            // isoV(0,2) = V(startIndex,2);

            int i = 0;
            for(auto element: pathPrint) {
                //cout << "element " << element << endl;
                //cout << "i " << i << endl;
                isoV(i,0) = V(element,0);
                isoV(i,1) = V(element,1);
                isoV(i,2) = V(element,2);
                i++;
            }
            int k = 0;
            for (int j = 0; j < pathPrint.size()-1; j++) {
                //cout << "j " << j << endl;
                //cout << "k " << k << endl;
                isoE(k,0) = j;
                isoE(k,1) = j+1;
                k++;
            }
            // cout << "Coloring2" << endl;
            //
            // C(startIndex, 0) = 0;
            // C(startIndex, 1) = 0;
            // C(startIndex, 2) = 0;
            //
            // for(auto element: pathPrint) {
            //
            //     C(element, 0) = 1;
            //     C(element, 1) = 0.4;
            //     C(element, 2) = 0.4;
            //
            // }
            // C(endIndex, 0) = 0;
            // C(endIndex, 1) = 0;
            // C(endIndex, 2) = 0;

            igl::opengl::glfw::Viewer viewer;

            viewer.data().set_edges(isoV,isoE,Eigen::RowVector3d(1,0,0));
            viewer.data().show_overlay_depth = false;

            viewer.data().set_mesh(V, F);
            //viewer.data().set_colors(C);
            viewer.launch();



        }
        else if (program == 3) {

            cout << "Applying Dijkstra (Needed further)" << endl;

            std::vector<std::vector<double>> adj_list;
            //ROW i is the number of the vertice, and the elements of that row indicates the neighbors.
            igl::adjacency_list(F,adj_list);

            std::vector<std::vector<float>> dijkstraMatrix;

            for(int i = 0; i < V.rows(); i++)
            {
                //cout << " Dijkstra number : " << i << "out of " << V.rows()  << endl;

                    //APPLY DIJKSTRA

                    std::vector<float> result = shortestPath(i, V.rows(), adj_list);

                    dijkstraMatrix.push_back(result);

            }

            //FIND FARTHEST

            std::vector<int> fpsSampleArray;

            int fartx = 0;
            int farty = 0;

            double maxDijkValue = 0.0;

            for(int i = 0; i < V.rows(); i++)
            {
                for(int j = 0; j < V.rows(); j++)
                {
                    if(dijkstraMatrix[i][j] > maxDijkValue) {

                        maxDijkValue = dijkstraMatrix[i][j];
                        fartx = i;
                        farty = j;

                    }

                }

            }

            fpsSampleArray.push_back(fartx);

            fpsSampleArray.push_back(farty);

            cout << "FPS SAMPLING STARTS" << endl;


            double farthestPointsDataDistance = 0;

            int farthestPointsDataNode = 0;


            for(int i = 0; i < 98; i++)
            {
                //cout << "sampling number: " << i << endl;
                float minValue = INF;//COMPARE BEFORE

                int minVertice = -1;

                int assocWith = -1;

                float maxValue = 0;//COMPARE AFTER

                int selectedVertice = -1;

                //TRAVERSE ALL VERTICES EXCEPT SAMPLED

                for(int j = 0; j < V.rows(); j++)
                {
                    if(find(fpsSampleArray.begin(), fpsSampleArray.end(), j) == fpsSampleArray.end()) { //IF NOT SAMPLED ALREADY



                        for(auto interator: fpsSampleArray) {
                            //cout << "-----" <<  endl;
                            //cout << "number j: " << j << endl;
                            //cout << "interator : " << interator << endl;

                            float temp = dijkstraMatrix[j][interator];

                            //cout << "temp j: " << temp << endl;
                            //cout << "minValue : " << minValue << endl;

                            if(temp < minValue) {
                                minValue = temp;
                                minVertice = j;
                                assocWith = interator;
                                //cout << "*new min vertice: " << minVertice << endl;
                            }
                            //cout << "-----" <<  endl;
                        }

                        //ASSOC WITH ONE OF THE SAMPLED, NOW HOLD BIGGEST ONE HERE

                        if(minValue > maxValue) {

                            //cout << "==new selected vertice: " << minVertice << endl;
                            //cout << "maxValue: " << maxValue << endl;
                            //cout << "minValue : " << minValue << endl;
                            selectedVertice = minVertice;

                            maxValue = minValue;

                            if(maxValue > farthestPointsDataDistance) {
                                farthestPointsDataDistance = maxValue;
                                farthestPointsDataNode = selectedVertice;
                                //farthestPointsDataInter = interator;
                            }
                        }

                        minValue = INF;//COMPARE BEFORE

                    }



                }

                //PUSH NEW SAMPLED INDEX
                fpsSampleArray.push_back(selectedVertice);


            }

            farthestPointsDataDistance = maxDijkValue;

            farthestPointsDataNode = fartx;

            // std::cout << "*****" << std::endl;

            // std::cout << "Max Dist: " << farthestPointsDataDistance << std::endl;

            // std::cout << "Max Node: " << farthestPointsDataNode << std::endl;

            // std::cout << "*****" << std::endl;

            std::cout << "FPS SAMPLING ENDED" << std::endl;






            //FPS DONE, ISO CURVE TIME

            //farthestPointsData holds 1: distance, 2: start point, 3: end point


            std::vector<nodeDistPair> sampledPairs;


            std::vector <nodeDistPair> sortedDijkForKIso;

            for(int i = 0; i < V.rows(); i++) {

                sortedDijkForKIso.push_back(make_pair(dijkstraMatrix[farthestPointsDataNode][i], i));

            }

        //RADII CALCULATION (SECTION 4.4)

            double delta = 0.02;

            double f_ls = 16;
            //std::cout << "fls: " << f_ls <<std::endl;

            double d_max = farthestPointsDataDistance;
            //std::cout << "dmax: " << d_max <<std::endl;


            double dfrac = d_max / f_ls;
            //std::cout << "dfrac: " << dfrac <<std::endl;



            std::vector<double> radiiList;

            for(int j = 0; j < 100; j++) {


                radiiList.push_back(((j+1)*delta)*((j+1)*delta)*((j+1)*delta)*((j+1)*delta) * dfrac);



            }



            std::cout << "LINE SEGMENT LENGTH COMPUTATION" << std::endl;
            //LINE SEGMENT LENGTH COMPUTATION
            //fpsSampleArray

            //FORLINE
            Eigen::VectorXd lineLenValues = Eigen::VectorXd::Zero(V.rows());

            double lineValues[100];

            for(int i = 0; i < 100; i++)
            {
                for(int j = 0; j < V.rows(); j++)
                {
                    if (F(j,0) == fpsSampleArray[i] || F(j,1) == fpsSampleArray[i] || F(j,2) == fpsSampleArray[i]) {
                        /* code */

                        //FOUND POINTS, F[j], APPLY ALGO

                        double r1 = *std::lower_bound(radiiList.begin(), radiiList.end(), dijkstraMatrix[farthestPointsDataNode][F(j,1)]);

                        double r2 = *std::lower_bound(radiiList.begin(), radiiList.end(), dijkstraMatrix[farthestPointsDataNode][F(j,2)]);

                        double alpha1 = std::abs(r1 - dijkstraMatrix[farthestPointsDataNode][F(j,0)]) / std::abs(dijkstraMatrix[farthestPointsDataNode][F(j,1)] - dijkstraMatrix[farthestPointsDataNode][F(j,0)]);

                        double alpha2 = std::abs(r2 - dijkstraMatrix[farthestPointsDataNode][F(j,0)]) / std::abs(dijkstraMatrix[farthestPointsDataNode][F(j,2)] - dijkstraMatrix[farthestPointsDataNode][F(j,0)]);

                        double p1x = (1 - alpha1) * V(F(j,0),0) + alpha1 * V(F(j,1),0);

                        double p1y = (1 - alpha1) * V(F(j,0),1) + alpha1 * V(F(j,1),1);

                        double p1z = (1 - alpha1) * V(F(j,0),2) + alpha1 * V(F(j,1),2);


                        double p2x = (1 - alpha2) * V(F(j,0),0) + alpha2 * V(F(j,2),0);

                        double p2y = (1 - alpha2) * V(F(j,0),1) + alpha2 * V(F(j,2),1);

                        double p2z = (1 - alpha2) * V(F(j,0),2) + alpha2 * V(F(j,2),2);

                        //Eigen::MatrixXd p1 = (1 - alpha1) * V[F(j,0)] + alpha1 * V[F(j,1)];

                        //Eigen::MatrixXd p2 = (1 - alpha2) * V[F(j,0)] + alpha2 * V[F(j,2)];

                        double distance = sqrt( (p1x-p2x)*(p1x-p2x) + (p1y-p2y)*(p1y-p2y) + (p1z-p2z)*(p1z-p2z) );

                        lineValues[i] = distance;
                        break;
                    }



                }

            }



            std::sort(sortedDijkForKIso.begin(), sortedDijkForKIso.end());

            Eigen::VectorXd isoValues = Eigen::VectorXd::Zero(V.rows());



            for(int i = 0; i < V.rows(); i++)
            {
                for(int j = 0; j < 100; j++)
                {
                    if (sortedDijkForKIso[i].first > radiiList[j]) {
                        /* code */
                        if (sortedDijkForKIso[i].first > farthestPointsDataDistance) {
                            /* code */
                            isoValues(sortedDijkForKIso[i].second,0) = farthestPointsDataDistance;
                            lineLenValues(sortedDijkForKIso[i].second, 0) = lineValues[farthestPointsDataNode];
                            break;
                        }


                        continue;
                    }
                    else
                    {
                        if (std::abs(sortedDijkForKIso[i].first - radiiList[j]) < std::abs(sortedDijkForKIso[i].first - radiiList[j-1]) ) {
                            isoValues(sortedDijkForKIso[i].second,0) = radiiList[j];
                            lineLenValues(sortedDijkForKIso[i].second, 0) = lineValues[j];
                            break;
                        }
                        else
                        {
                            /* code */
                            isoValues(sortedDijkForKIso[i].second,0) = radiiList[j-1];
                            lineLenValues(sortedDijkForKIso[i].second, 0) = lineValues[j-1];
                            break;
                        }


                    }


                }

            }

            //std::cout << "----lineLen" << std::endl;

            //std::cout << lineLenValues << std::endl;

            //std::cout << "----lineLen" << std::endl;

            Eigen::MatrixXd isoV;
            Eigen::MatrixXi isoE;

            const Eigen::VectorXd* z;
            const Eigen::VectorXd* zColor;

            z = &isoValues;
            zColor = &lineLenValues;

            igl::isolines(V, F, *z, 100, isoV, isoE);


            // Plot the mesh
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(V, F);

            viewer.data().set_edges(isoV,isoE,Eigen::RowVector3d(0,0,0));
            Eigen::MatrixXd colors;
            igl::jet(*zColor, true, colors);
            viewer.data().set_colors(colors);
            std::cout << "$$$$$$$$$$" << std::endl;
            std::cout << "IMPORTANT: To check isolines better, press 'l' to disable wireframe and 'f' to enable face based modeling." << std::endl;
            std::cout << "$$$$$$$$$$" << std::endl;
            std::cout << "$$$$$$$$$$" << std::endl;
            std::cout << "IMPORTANT: To close the viewer, press ESC while using it." << std::endl;
            std::cout << "$$$$$$$$$$" << std::endl;

            viewer.launch();




        }
        else if (program == 4) {
            std::cout << "Take 2 seed vertices as input and draw the shortest path 'p' between them. Then color the remainingvertices based on their distances to 'p'." << std::endl;


            int startIndex;
            int endIndex;

            std::cout << "Please enter the index # of the starting vertice out of " << V.rows() << " vertex" << std::endl;
            cin >> startIndex;
            std::cout << "Please enter the index # of the ending vertice out of " << V.rows() << " vertex" << std::endl;
            cin >> endIndex;
            //TODO:HIGHLIGHT PATH, ALSO PRINT PATH

            Eigen::MatrixXd C; //COLOR MATRIX

            C = Eigen::MatrixXd::Constant(V.rows(),3,1); //WHITE ALL



            Eigen::Matrix<float,Eigen::Dynamic,1> T;
            Eigen::Matrix<double,Eigen::Dynamic,1> sD;

            igl::project_to_line(V,V.row(startIndex).eval(),V.row(endIndex).eval(),T,sD);

            //std::cout << sD << std::endl;



            igl::opengl::glfw::Viewer viewer;


            //draw line


            C = Eigen::MatrixXd::Constant(V.rows(),3,1); //WHITE ALL

            cout << "Applying Dijkstra" << endl;
            std::vector<std::vector<double>> adj_list;
            //ROW i is the number of the vertice, and the elements of that row indicates the neighbors.
            igl::adjacency_list(F,adj_list);

            std::vector<int >path = shortestPathPrint(startIndex, endIndex, V.rows(), adj_list);

        //    std::cout << "####" << std::endl;
        //    for (int i=0;i < V.rows();i++)
        //        std::cout << path[i] << ' ';
        //    std::cout << "####" << std::endl;


            cout << "Coloring" << endl;
            std::vector<int> pathPrint;
            pathPrint.push_back(startIndex);
            printPath(path, endIndex, pathPrint);


            Eigen::MatrixXd isoV = Eigen::MatrixXd::Zero(pathPrint.size(),3);
            Eigen::MatrixXi isoE = Eigen::MatrixXi::Zero(pathPrint.size()-1,2);

            //cout << "size " << pathPrint.size()+1 << endl;

            // std::cout << "####" << std::endl;
            // for (auto iter: pathPrint)
            //     std::cout << iter << ' ';
            // std::cout << "####" << std::endl;

            // isoV(0,0) = V(startIndex,0);
            // isoV(0,1) = V(startIndex,1);
            // isoV(0,2) = V(startIndex,2);

            //Eigen::VectorXi FS,VT,FT;

            //Eigen::Matrix<int,Eigen::Dynamic,1> VS;

            //Eigen::MatrixXd VS = Eigen::MatrixXd::Zero(pathPrint.size(),1);

            int i = 0;
            for(auto element: pathPrint) {
                //cout << "element " << element << endl;
                //cout << "i " << i << endl;
                isoV(i,0) = V(element,0);
                isoV(i,1) = V(element,1);
                isoV(i,2) = V(element,2);

                //VS(i)=element;

                i++;
            }
            int k = 0;
            for (int j = 0; j < pathPrint.size()-1; j++) {
                //cout << "j " << j << endl;
                //cout << "k " << k << endl;
                isoE(k,0) = j;
                isoE(k,1) = j+1;
                k++;
            }
            // cout << "Coloring2" << endl;
            //
            // C(startIndex, 0) = 0;
            // C(startIndex, 1) = 0;
            // C(startIndex, 2) = 0;
            //
            // for(auto element: pathPrint) {
            //
            //     C(element, 0) = 1;
            //     C(element, 1) = 0.4;
            //     C(element, 2) = 0.4;
            //
            // }
            // C(endIndex, 0) = 0;
            // C(endIndex, 1) = 0;
            // C(endIndex, 2) = 0;


            viewer.data().set_edges(isoV,isoE,Eigen::RowVector3d(1,0,0));
            viewer.data().show_overlay_depth = false;

            

            //Eigen::VectorXd d;

            //VT.setLinSpaced(V.rows(),0,V.rows()-1);

            //igl::exact_geodesic(V,F,VS,FS,VT,FT,d);

            //cout << "??????? " <<  endl;



            viewer.data().set_mesh(V, F);
            viewer.data().set_colors(sD);
            viewer.launch();



        }
        else if (program == 5) {
            std::cout << "See you" << std::endl;
            break;
        }
    }
    return 0;
}
