CC=g++
BOOST_ROOT?=/usr/local
RDKIT_ROOT?=/usr/local
MONGO_CXX_DRIVER_ROOT?=/usr/local

bin/lbvs: obj/main.o obj/compound_database.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${MONGO_CXX_DRIVER_ROOT}/lib64 -lmongocxx -lbsoncxx -L${RDKIT_ROOT}/lib -lRDKitDescriptors -lRDKitFingerprints -lRDKitFileParsers -lRDKitDepictor -lRDKitMolTransforms -lRDKitSubstructMatch -lRDKitSmilesParse -lRDKitAlignment -lRDKitGraphMol -lRDKitRDGeometryLib -lRDKitRDGeneral -L${BOOST_ROOT}/lib -lboost_date_time

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -Wno-deprecated-declarations -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit -I${MONGO_CXX_DRIVER_ROOT}/include/mongocxx/v_noabi -I${MONGO_CXX_DRIVER_ROOT}/include/bsoncxx/v_noabi

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${BOOST_ROOT}/include

clean:
	rm -f bin/lbvs obj/*.o
