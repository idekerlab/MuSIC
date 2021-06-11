#include <iostream>
#include <stdlib.h>
#include "dagConstruct.h"


int main(int argc, char* argv[]) {

  if (argc < 3) {
    cout << "Needs 2 arguments: 1) file with clusters of terminal nodes (i.e. genes) listed, one cluster per line, comma separated" << endl;
    cout << "2) number of terminal nodes in DAG (input must be at least as large as number of terminal nodes in DAG)" << endl;
    cout << "3) optional - number of clusters in file" << endl;
    cout << "4) Category name of terminal nodes (i.e. gene, patient. etc.). Optional with default = gene" << endl;
    return 0;
  }

  map<string, unsigned> nodeNamesToIDs;
  unsigned nextID = 0;
  vector<validClusterBitset> validClusters;
  string clustersFile = argv[1];
  unsigned numTerminalNodes = stoi(argv[2]);
  unsigned nextClustID = numTerminalNodes;
  if (argc >= 4) {
    validClusters.reserve(stoi(argv[3]));
  }
  DAGraph ontology;
  if (argc >= 5) {
    ontology.setTerminalName(argv[4]);
  }

  time_t start, end;
  time(&start);
  //cout << "Loading clusters" << endl;
  string line;
  ifstream file(clustersFile.c_str());
  if (file.is_open()) {
    while (file.good()) {
      getline(file,line);
      boost::dynamic_bitset<unsigned long> cluster(numTerminalNodes);
      vector<string> tokens;
      Utils::Tokenize(line, tokens, ",");
      for (unsigned i = 0; i < tokens.size(); ++i) {
	map<string,unsigned int>::iterator nodeIDIt = nodeNamesToIDs.find(tokens[i]);
	if (nodeIDIt == nodeNamesToIDs.end()) {
	  cluster[nextID] = 1;
	  nodeNamesToIDs.insert(make_pair(tokens[i], nextID));
	  ontology.addNode(tokens[i], nodeNamesToIDs);
	  ++nextID;
	  if (nextID > numTerminalNodes) {
	    cout << "Number of terminal nodes greater than input - exiting" << endl;
	  }
	} else {
	  cluster[nodeIDIt->second] = 1;
	}
      }
      if (tokens.size() > 0) {
	validClusters.push_back(validClusterBitset(cluster, nextClustID, (static_cast<double>(tokens.size()) / numTerminalNodes) / 2));
	++nextClustID;
      }
    }
  }
  
  //cout << "Actual number of terminal nodes is: " << nodeNamesToIDs.size() << endl;
  time(&end);
  double dif = difftime(end,start);
  cout << "loading clusters took " << dif << " seconds" << endl;
  
  // Create ontology from other networks
  //cout << "Creating DAG" << endl;
  time (&start);
  dagConstruct::clustersToDAG(validClusters, ontology, numTerminalNodes);
  time (&end);
  dif = difftime(end,start);
  cout << "Creating DAG took " << dif << " seconds" << endl;

  //cout << "Ontology is: " << endl;
  for(map< pair<unsigned,unsigned>, string >::iterator edgesIt = ontology.edgesBegin(); edgesIt != ontology.edgesEnd(); ++edgesIt) {
    cout << ontology.getName(edgesIt->first.first) << "\t" << ontology.getName(edgesIt->first.second) << "\t" << edgesIt->second << "\t" << endl;
  }
  return 1;
}
