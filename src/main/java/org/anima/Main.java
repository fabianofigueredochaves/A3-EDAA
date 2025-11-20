package org.anima;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Main {

    // Estrutura para guardar o resultado de um algoritmo
    static class Resultado {
        List<Integer> caminho;
        double custo;
        int nosExpandidos;
        double tempoMs;

        Resultado(List<Integer> caminho, double custo, int nosExpandidos, double tempoMs) {
            this.caminho = caminho;
            this.custo = custo;
            this.nosExpandidos = nosExpandidos;
            this.tempoMs = tempoMs;
        }
    }

    // Nó para filas de prioridade (Dijkstra, Greedy, A*)
    static class NoPQ {
        int v;
        double prioridade; // pode ser dist, f = g + h, ou apenas h

        NoPQ(int v, double prioridade) {
            this.v = v;
            this.prioridade = prioridade;
        }
    }

    // ---------- Leitura da matriz de adjacência ----------
    private static int[][] lerMatrizAdjacencia(String arquivo) throws IOException {
        List<String> linhas = Files.readAllLines(Paths.get(arquivo));
        List<int[]> linhasValores = new ArrayList<>();

        for (String linha : linhas) {
            linha = linha.trim();
            if (linha.isEmpty()) continue;

            String[] partes = linha.split("\\s+");
            int[] valores = new int[partes.length];
            for (int i = 0; i < partes.length; i++) {
                valores[i] = Integer.parseInt(partes[i]);
            }
            linhasValores.add(valores);
        }
        int n = linhasValores.size();
        int[][] matriz = new int[n][n];
        for (int i = 0; i < n; i++) {
            int[] linhaVals = linhasValores.get(i);
            if (linhaVals.length != n) {
                throw new IOException("Matriz de adjacência não é quadrada. Linha " + i);
            }
            System.arraycopy(linhaVals, 0, matriz[i], 0, n);
        }
        return matriz;
    }

    // ---------- Conversão coordenada (linha,coluna) <-> índice ----------
    private static int coordParaIndice(String coordStr, int cols) {
        String[] partes = coordStr.split(",");
        int linha = Integer.parseInt(partes[0].trim());
        int coluna = Integer.parseInt(partes[1].trim());
        return linha * cols + coluna;
    }

    private static String indiceParaCoord(int idx, int cols) {
        int linha = idx / cols;
        int coluna = idx % cols;
        return "(" + linha + "," + coluna + ")";
    }

    // ---------- Funções auxiliares de caminho, custo, saída ----------
    private static List<Integer> reconstruirCaminho(int[] pai, int origem, int destino) {
        if (origem == destino && pai[destino] == -1) {
            // Caminho trivial
            List<Integer> caminho = new ArrayList<>();
            caminho.add(origem);
            return caminho;
        }
        if (pai[destino] == -1) {
            return null; // sem caminho
        }

        List<Integer> caminho = new ArrayList<>();
        int atual = destino;
        while (atual != -1) {
            caminho.add(atual);
            atual = pai[atual];
        }
        Collections.reverse(caminho);
        return caminho;
    }

    private static double calcularCusto(List<Integer> caminho, int[][] grafo) {
        if (caminho == null || caminho.size() < 2) return 0.0;
        double custo = 0.0;
        for (int i = 0; i < caminho.size() - 1; i++) {
            int u = caminho.get(i);
            int v = caminho.get(i + 1);
            custo += grafo[u][v];
        }
        return custo;
    }

    private static String formatarCaminho(List<Integer> caminho, int cols) {
        if (caminho == null || caminho.isEmpty()) return "";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < caminho.size(); i++) {
            sb.append(indiceParaCoord(caminho.get(i), cols));
            if (i < caminho.size() - 1) {
                sb.append(" -> ");
            }
        }
        return sb.toString();
    }

    private static void escreverResultado(String arquivoSaida,
                                          String nomeAlgoritmo,
                                          String heuristica,
                                          String origemStr,
                                          String destinoStr,
                                          Resultado res,
                                          int cols) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(arquivoSaida))) {
            pw.println("ALGORITIMO: " + nomeAlgoritmo);
            pw.println("HEURISTICA: " + (heuristica == null ? "" : heuristica));
            pw.println("ORIGEM: (" + origemStr + ")");
            pw.println("DESTINO: (" + destinoStr + ")");
            String caminhoStr = (res != null && res.caminho != null)
                    ? formatarCaminho(res.caminho, cols)
                    : "";
            pw.println("CAMINHO: " + caminhoStr);
            String custoStr = (res != null && res.caminho != null)
                    ? String.valueOf(res.custo)
                    : "";
            String nosStr = (res != null) ? String.valueOf(res.nosExpandidos) : "";
            String tempoStr = (res != null)
                    ? String.format(Locale.US, "%.2f", res.tempoMs)
                    : "";

            pw.println("CUSTO: " + custoStr);
            pw.println("NOS EXPANDIDOS: " + nosStr);
            pw.println("TEMPO (ms): " + tempoStr);
        }
    }

    // ---------- BFS (Breadth-First Search) ----------
    public static Resultado bfs(int[][] grafo, int origem, int destino) {
        int n = grafo.length;
        boolean[] visitado = new boolean[n];
        int[] pai = new int[n];
        Arrays.fill(pai, -1);
        Queue<Integer> fila = new ArrayDeque<>();
        visitado[origem] = true;
        fila.add(origem);
        int nosExpandidos = 0;
        long inicio = System.nanoTime();
        while (!fila.isEmpty()) {
            int u = fila.poll();
            nosExpandidos++;
            if (u == destino) break;

            for (int v = 0; v < n; v++) {
                if (grafo[u][v] > 0 && !visitado[v]) {
                    visitado[v] = true;
                    pai[v] = u;
                    fila.add(v);
                }
            }
        }
        long fim = System.nanoTime();
        double tempoMs = (fim - inicio) / 1_000_000.0;
        List<Integer> caminho = (destino == origem || pai[destino] != -1)
                ? reconstruirCaminho(pai, origem, destino)
                : null;

        double custo = calcularCusto(caminho, grafo);
        return new Resultado(caminho, custo, nosExpandidos, tempoMs);

    }

    // ---------- DFS (Depth-First Search) iterativo ----------
    public static Resultado dfs(int[][] grafo, int origem, int destino) {
        int n = grafo.length;
        boolean[] visitado = new boolean[n];
        int[] pai = new int[n];
        Arrays.fill(pai, -1);
        Deque<Integer> pilha = new ArrayDeque<>();
        pilha.push(origem);
        int nosExpandidos = 0;
        long inicio = System.nanoTime();
        while (!pilha.isEmpty()) {
            int u = pilha.pop();
            if (visitado[u]) continue;
            visitado[u] = true;
            nosExpandidos++;
            if (u == destino) break;

            // Pode influenciar o caminho; aqui vai de 0 até n-1
            for (int v = 0; v < n; v++) {
                if (grafo[u][v] > 0 && !visitado[v]) {
                    if (pai[v] == -1) pai[v] = u;
                    pilha.push(v);
                }
            }
        }
        long fim = System.nanoTime();
        double tempoMs = (fim - inicio) / 1_000_000.0;
        List<Integer> caminho = (destino == origem || pai[destino] != -1)
                ? reconstruirCaminho(pai, origem, destino)
                : null;

        double custo = calcularCusto(caminho, grafo);
        return new Resultado(caminho, custo, nosExpandidos, tempoMs);
    }

    // ---------- Dijkstra ----------
    public static Resultado dijkstra(int[][] grafo, int origem, int destino) {
        int n = grafo.length;
        double[] dist = new double[n];
        int[] pai = new int[n];
        boolean[] visitado = new boolean[n];
        Arrays.fill(dist, Double.POSITIVE_INFINITY);
        Arrays.fill(pai, -1);
        dist[origem] = 0.0;
        PriorityQueue<NoPQ> fila = new PriorityQueue<>(Comparator.comparingDouble(a -> a.prioridade));
        fila.add(new NoPQ(origem, 0.0));
        int nosExpandidos = 0;
        long inicio = System.nanoTime();
        while (!fila.isEmpty()) {
            NoPQ atual = fila.poll();
            int u = atual.v;
            if (visitado[u]) continue;
            visitado[u] = true;
            nosExpandidos++;
            if (u == destino) break;

            for (int v = 0; v < n; v++) {
                if (grafo[u][v] > 0 && !visitado[v]) {
                    double custoAresta = grafo[u][v];
                    double novaDist = dist[u] + custoAresta;
                    if (novaDist < dist[v]) {
                        dist[v] = novaDist;
                        pai[v] = u;
                        fila.add(new NoPQ(v, novaDist));
                    }
                }
            }
        }
        long fim = System.nanoTime();
        double tempoMs = (fim - inicio) / 1_000_000.0;
        if (Double.isInfinite(dist[destino])) {
            return new Resultado(null, 0.0, nosExpandidos, tempoMs);
        }

        List<Integer> caminho = reconstruirCaminho(pai, origem, destino);
        double custo = dist[destino];
        return new Resultado(caminho, custo, nosExpandidos, tempoMs);
    }

    // ---------- Heurísticas ----------
    private static double heuristicaManhattan(int v, int destino, int cols) {
        int vx = v / cols;
        int vy = v % cols;
        int dx = destino / cols;
        int dy = destino % cols;
        return Math.abs(vx - dx) + Math.abs(vy - dy);
    }

    private static double heuristicaEuclidiana(int v, int destino, int cols) {
        int vx = v / cols;
        int vy = v % cols;
        int dx = destino / cols;
        int dy = destino % cols;
        int difx = vx - dx;
        int dify = vy - dy;
        return Math.sqrt(difx * difx + dify * dify);
    }

    // ---------- Greedy Best-First Search ----------
    public static Resultado greedyBestFirst(int[][] grafo, int origem, int destino,
                                            int cols, boolean manhattan) {
        int n = grafo.length;
        boolean[] visitado = new boolean[n];
        int[] pai = new int[n];
        Arrays.fill(pai, -1);
        PriorityQueue<NoPQ> fila = new PriorityQueue<>(Comparator.comparingDouble(a -> a.prioridade));
        double hOrigem = manhattan
                ? heuristicaManhattan(origem, destino, cols)
                : heuristicaEuclidiana(origem, destino, cols);
        fila.add(new NoPQ(origem, hOrigem));
        int nosExpandidos = 0;
        long inicio = System.nanoTime();
        while (!fila.isEmpty()) {
            NoPQ atual = fila.poll();
            int u = atual.v;
            if (visitado[u]) continue;
            visitado[u] = true;
            nosExpandidos++;
            if (u == destino) break;

            for (int v = 0; v < n; v++) {
                if (grafo[u][v] > 0 && !visitado[v]) {
                    if (pai[v] == -1) pai[v] = u;
                    double h = manhattan
                            ? heuristicaManhattan(v, destino, cols)
                            : heuristicaEuclidiana(v, destino, cols);
                    fila.add(new NoPQ(v, h)); // só a heurística manda
                }
            }
        }
        long fim = System.nanoTime();
        double tempoMs = (fim - inicio) / 1_000_000.0;

        List<Integer> caminho = (origem == destino || pai[destino] != -1)
                ? reconstruirCaminho(pai, origem, destino)
                : null;
        double custo = calcularCusto(caminho, grafo);
        return new Resultado(caminho, custo, nosExpandidos, tempoMs);
    }

    // ---------- A* (A-estrela) ----------
    public static Resultado aStar(int[][] grafo, int origem, int destino,
                                  int cols, boolean manhattan) {
        int n = grafo.length;
        double[] g = new double[n];
        double[] f = new double[n];
        int[] pai = new int[n];
        boolean[] fechado = new boolean[n];
        Arrays.fill(g, Double.POSITIVE_INFINITY);
        Arrays.fill(f, Double.POSITIVE_INFINITY);
        Arrays.fill(pai, -1);
        g[origem] = 0.0;
        double hOrigem = manhattan
                ? Main.heuristicaManhattan(origem, destino, cols)
                : heuristicaEuclidiana(origem, destino, cols);
        f[origem] = hOrigem;
        PriorityQueue<NoPQ> fila = new PriorityQueue<>(Comparator.comparingDouble(a -> a.prioridade));
        fila.add(new NoPQ(origem, f[origem]));
        int nosExpandidos = 0;
        long inicio = System.nanoTime();
        while (!fila.isEmpty()) {
            NoPQ atual = fila.poll();
            int u = atual.v;
            if (fechado[u]) continue;
            fechado[u] = true;
            nosExpandidos++;

            if (u == destino) break;

            for (int v = 0; v < n; v++) {
                if (grafo[u][v] > 0 && !fechado[v]) {
                    double custoAresta = grafo[u][v];
                    double gTentativa = g[u] + custoAresta;
                    if (gTentativa < g[v]) {
                        pai[v] = u;
                        g[v] = gTentativa;
                        double h = manhattan
                                ? heuristicaManhattan(v, destino, cols)
                                : heuristicaEuclidiana(v, destino, cols);
                        f[v] = g[v] + h;
                        fila.add(new NoPQ(v, f[v]));
                    }
                }
            }
        }
        long fim = System.nanoTime();
        double tempoMs = (fim - inicio) / 1_000_000.0;
        if (Double.isInfinite(g[destino])) {
            return new Resultado(null, 0.0, nosExpandidos, tempoMs);
        }
        List<Integer> caminho = reconstruirCaminho(pai, origem, destino);
        double custo = g[destino];
        return new Resultado(caminho, custo, nosExpandidos, tempoMs);
    }

    public static void main(String[] args) {
        // 1. Definir diretórios de entrada e saída
        String dirEntrada = "matrizes";
        String dirResultado = "resultado";

        // 2. Lista dos arquivos de matriz para processar
        List<String> nomesArquivos = Arrays.asList(
                "matrix_4x4.txt",
                "matrix_16x16.txt",
                "matrix_32x32.txt",
                "matrix_64x64.txt"
        );

        // Garante que o diretório de resultados exista
        try {
            Files.createDirectories(Paths.get(dirResultado));
        } catch (IOException e) {
            System.err.println("Erro ao criar o diretório de resultado: " + e.getMessage());
            return; // Encerra se não conseguir criar a pasta
        }

        // 3. Loop para processar cada arquivo
        for (String nomeArquivo : nomesArquivos) {
            Path caminhoEntrada = Paths.get(dirEntrada, nomeArquivo);
            // Extrai o nome do arquivo sem a extensão .txt para usar na saída
            String baseNomeSaida = nomeArquivo.replace(".txt", "");
            Path baseCaminhoSaida = Paths.get(dirResultado, baseNomeSaida);

            System.out.println("----------------------------------------------------");
            System.out.println("Processando matriz: " + caminhoEntrada);

            try {
                int[][] grafo = Main.lerMatrizAdjacencia(caminhoEntrada.toString());
                int n = grafo.length;

                int cols = (int) Math.round(Math.sqrt(n));
                if (cols * cols != n) {
                    System.err.println("Aviso: número de vértices não é quadrado perfeito para " + nomeArquivo);
                    continue; // Pula para o próximo arquivo
                }

                // 4. Definir origem e destino padrão (canto superior esquerdo ao inferior direito)
                String origemStr = "0,0";
                String destinoStr = (cols - 1) + "," + (cols - 1);
                int origem = coordParaIndice(origemStr, cols);
                int destino = coordParaIndice(destinoStr, cols);

                System.out.println("Origem: (" + origemStr + ") -> Destino: (" + destinoStr + ")");


                // 5. Executa todos os algoritmos para a matriz atual
                Resultado resBfs = bfs(grafo, origem, destino);
                Resultado resDfs = dfs(grafo, origem, destino);
                Resultado resDijkstra = dijkstra(grafo, origem, destino);
                Resultado resGreedyManhattan = greedyBestFirst(grafo, origem, destino, cols, true);
                Resultado resGreedyEuclidiana = greedyBestFirst(grafo, origem, destino, cols, false);
                Resultado resAStarManhattan = aStar(grafo, origem, destino, cols, true);
                Resultado resAStarEuclidiana = aStar(grafo, origem, destino, cols, false);

                // 6. Grava os arquivos de saída no diretório 'resultado'
                escreverResultado(baseCaminhoSaida + ".bfs", "BFS", "", origemStr, destinoStr, resBfs, cols);
                escreverResultado(baseCaminhoSaida + ".dfs", "DFS", "", origemStr, destinoStr, resDfs, cols);
                escreverResultado(baseCaminhoSaida + ".dijkstra", "DIJKSTRA", "", origemStr, destinoStr, resDijkstra, cols);
                escreverResultado(baseCaminhoSaida + ".gbs.manhattan", "GREEDY BEST-FIRST-SEARCH", "Manhattan", origemStr, destinoStr, resGreedyManhattan, cols);
                escreverResultado(baseCaminhoSaida + ".gbs.euclidiana", "GREEDY BEST-FIRST-SEARCH", "Euclidiana", origemStr, destinoStr, resGreedyEuclidiana, cols);
                escreverResultado(baseCaminhoSaida + ".a.manhattan", "A*", "Manhattan", origemStr, destinoStr, resAStarManhattan, cols);
                escreverResultado(baseCaminhoSaida + ".a.euclidiana", "A*", "Euclidiana", origemStr, destinoStr, resAStarEuclidiana, cols);

                System.out.println("Execução para " + nomeArquivo + " concluída com sucesso.");

            } catch (Exception ex) {
                System.err.println("Erro ao processar o arquivo " + nomeArquivo + ": " + ex.getMessage());
                ex.printStackTrace();
            }
        }
        System.out.println("----------------------------------------------------");
        System.out.println("Todos os arquivos foram processados. Verifique a pasta '" + dirResultado + "'.");
    }
}
/*
    public static void main(String[] args) {

        if (args.length < 3) {
             System.err.println("Uso: java BuscaRotasGrafos <arquivo_entrada>, linha, coluna");
            System.exit(1);
        }
        String arquivoEntrada = args[0];
        String origemStr = args[1].replace("\"", "");
        String destinoStr = args[2].replace("\"", "");
        //String arquivoEntrada = "0 1 0 1 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 1 0 1 0";
        //String origemStr = "0,0";
        //String destinoStr ="2,2";

        try {
            int[][] grafo = Main.lerMatrizAdjacencia(arquivoEntrada);
            int n = grafo.length;

            // Assume grade N x N, logo raiz quadrada perfeita
            int cols = (int) Math.round(Math.sqrt(n));
            if (cols * cols != n) {
                System.err.println("Aviso: número de vértices não é quadrado perfeito. Verifique o arquivo.");
            }
            int origem = coordParaIndice(origemStr, cols);
            int destino = coordParaIndice(destinoStr, cols);

            // Executa algoritmos
            Resultado resBfs = bfs(grafo, origem, destino);
            Resultado resDfs = dfs(grafo, origem, destino);
            Resultado resDijkstra = dijkstra(grafo, origem, destino);
            Resultado resGreedyManhattan = greedyBestFirst(grafo, origem, destino, cols, true);
            Resultado resGreedyEuclidiana = greedyBestFirst(grafo, origem, destino, cols, false);
            Resultado resAStarManhattan = aStar(grafo, origem, destino, cols, true);
            Resultado resAStarEuclidiana = aStar(grafo, origem, destino, cols, false);

            // Grava arquivos de saída
            escreverResultado(arquivoEntrada + ".bfs", "BFS", "", origemStr, destinoStr, resBfs, cols);
            escreverResultado(arquivoEntrada + ".dfs", "DFS", "", origemStr, destinoStr, resDfs, cols);
            escreverResultado(arquivoEntrada + ".dijkstra", "DIJKSTRA", "", origemStr, destinoStr, resDijkstra, cols);
            escreverResultado(arquivoEntrada + ".gbs.manhattan", "GREEDY BEST-FIRST-SEARCH",
                    "Manhattan", origemStr, destinoStr, resGreedyManhattan, cols);
            escreverResultado(arquivoEntrada + ".gbs.euclidiana", "GREEDY BEST-FIRST-SEARCH",
                    "Euclidiana", origemStr, destinoStr, resGreedyEuclidiana, cols);
            escreverResultado(arquivoEntrada + ".a.manhattan", "A*", "Manhattan",
                             origemStr, destinoStr, resAStarManhattan, cols);
            escreverResultado(arquivoEntrada + ".a.euclidiana", "A*", "Euclidiana",
                             origemStr, destinoStr, resAStarEuclidiana, cols);
            System.out.println("Execução concluída. Arquivos de saída gerados ao lado do arquivo de entrada.");
        } catch (Exception ex) {
      //      throw new RuntimeException(ex);
        }
        //} catch (IOException e;) {
        // System.err.println("Erro de IO: " + e.getMessage());
    } //catch (Exception e) {
    //   System.err.println("Erro: " + e.getMessage());
    //   e.printStackTrace();
    //}
    //}

}




public class Main {



    // [FIM DA MODIFICAÇÃO]

    // ... (Cole aqui o restante das suas classes e métodos: Resultado, NoPQ, lerMatrizAdjacencia, etc.)
}
 */