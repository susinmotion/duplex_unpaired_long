secondtrie:  secondtrie.cpp node.cpp trie.cpp variants.cpp initialize.cpp leafdata.cpp 
	g++ -pg -g secondtrie.cpp node.cpp trie.cpp variants.cpp initialize.cpp leafdata.cpp -o secondtrie

clean: 
	  $(RM) secondtrie
