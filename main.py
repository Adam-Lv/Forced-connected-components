from graph import GrapheFortementConnexe as gfc, TestGraph


def main():
    test = TestGraph()
    for i, G in enumerate(test):
        print("------------------------------------------------------------------------------------------")
        print("")
        print("Graphe", i + 1, ":")
        test.demo(G)
        print("")

    # G = gfc(test[1])
    # G.kosaraju()
    # G.connexion()
    # G.dessiner_le_processus(True)


if __name__ == '__main__':
    main()
