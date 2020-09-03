# orthogonal-art-gallery

Na repozitorijumu se nalaze programi oag.cpp i greedy.cpp koji za problem ortogonalne galerije (engl. Orthgonal Art Gallery) pronalaze stražare koji pokrivaju poligon koji se da programu kao ulaz. Program koristi heuristički algoritam opisan u master radu "Kompleksost problem ortogonalne galserije".

Biblioteke potrebne za rad porgrama oag.cpp su:
- cgal (verzija 5.0.2)
- minisat (verzija 2.2)

Verzija c++ kompajlera: gcc (Ubuntu 9.3.0-10ubuntu2) 9.3.0

Komanda kojom se pokreće oag.cpp:
g++ -g oag.cpp -o oag.out -lgmp -lmpfr -lminisat && clear && ./oag.out
Komanda kojom se pokreće greedy.cpp:
g++ -g greedy.cpp -o greedy.out -lgmp -lmpfr && clear && ./greedy.out

Takođe, na repozitorijumu se nalaze rezultati dobijeni nad instancama iz rada "An Exact and Efficient Algorithm for the Orthogonal Art Gallery Problem".
