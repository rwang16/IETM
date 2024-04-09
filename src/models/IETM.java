package models;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import utility.FuncUtils;
public class IETM
{
    public double alpha; // Hyper-parameter alpha
    public double beta; // Hyper-parameter beta
    public double weights; // Hyper-parameter weights
    public int K; // Number of topics

    public int topWords; // Number of most probable words for each topic

    public List<List<Integer>> ShortCorpus; // Word ID-based corpus
    public List<List<Integer>> pseudoCorpus;

    // in the corpus
    public int numShortDocuments; // Number of documents in the corpus
    public int numWordsInShortCorpus; // Number of words in the corpus

    // in the corpus
    public int numPseudoDocuments; // Number of documents in the corpus
    public int numWordsInPseudoCorpus; // Number of words in the corpus

    public HashMap<String, Integer> word2IdVocabulary; // Vocabulary to get ID
    // given a word
    public HashMap<Integer, String> id2WordVocabulary; // Vocabulary to get word
    // given an ID
    public int V; // The number of word types in the corpus

    /**
     * topic assignments for each word.
     */
    int[][] Pz;

    int[][] sz;
    /**
     * nw[i][j] number of instances of word i (term?) assigned to topic j.
     */
    int[][] nw;

    /**
     * ngd[i][j] number of words in pseudo document i assigned to topic j.
     */
    double[][] npd;
    /**
     * nsd[i][j] number of words in short document i assigned to topic j.
     */
    int[][] nsd;

    /**
     * nwsum[j] total number of words assigned to topic j.
     */
    int[] nwsum;

    /**
     * npdsum[i] total number of words in pseudo document i.
     */
    int[] npdsum;

    /**
     * nssum[i] total number of words in document i.
     */
    int[] nsdsum;

    /**
     * cumulative statistics of theta
     */
    double[][] thetasum;


    /**
     * cumulative statistics of phi
     */
    double[][] phisum;

    double[][] phi;
    double[] pz;
    double[][] pdz;

    /**
     * size of statistics
     */
    int numstats;
    /**
     * sampling lag (?)
     */
    private static int THIN_INTERVAL = 20;

    /**
     * burn-in period
     */
    private static int BURN_IN = 100;

    /**
     * max iterations
     */
    private static int ITERATIONS = 1000;

    private static int dispcol = 0;
    /**
     * sample lag (if -1 only one sample taken)
     */
    private static int SAMPLE_LAG;
    // Double array used to sample a topic
    public double[] multiPros;

    // Path to the directory containing the corpus
    public String folderPath;
    // Path to the topic modeling corpus
    public String ShortCorpusPath;
    public String pseudoCorpusPath;

    public String expName = "IETM";
    public String orgExpName = "IETM";


    public double initTime = 0;
    public double iterTime = 0;

    public IETM(String pathToShortCorpus, String pathToPseudoCorpus, String pathToresult, int inNumTopics,
                double inAlpha, double inBeta, int inNumIterations, int inTopWords)
            throws Exception
    {
        this(pathToShortCorpus, pathToPseudoCorpus,
                pathToresult, inNumTopics, inAlpha, inBeta, inNumIterations,
                inTopWords, "IETM");
    }

    public IETM(String pathToShortCorpus, String pathToPseudoCorpus, String pathToresult, int inNumTopics,
                double inAlpha, double inBeta, int inNumIterations, int inTopWords,
                String inExpName)
            throws Exception
    {

        alpha = inAlpha;
        beta = inBeta;
        K = inNumTopics;
        ITERATIONS = inNumIterations;
        topWords = inTopWords;

        expName = inExpName;
        orgExpName = expName;
        ShortCorpusPath = pathToShortCorpus;
        pseudoCorpusPath = pathToPseudoCorpus;

        folderPath = pathToresult;
        File dir = new File(folderPath);
        if (!dir.exists()) {
            dir.mkdirs();
        }
        System.out.println("Reading topic modeling corpus: " + ShortCorpusPath);
        System.out.println("Reading pseudo corpus: " + pseudoCorpusPath);

        word2IdVocabulary = new HashMap<String, Integer>();
        id2WordVocabulary = new HashMap<Integer, String>();
        ShortCorpus = new ArrayList<List<Integer>>();
        pseudoCorpus = new ArrayList<List<Integer>>();
        numShortDocuments = 0;
        numWordsInShortCorpus = 0;
        numPseudoDocuments = 0;
        numWordsInPseudoCorpus = 0;

        BufferedReader br = null;
        int indexWord = -1;
        try {
            br = new BufferedReader(new FileReader(ShortCorpusPath));
            for (String doc; (doc = br.readLine()) != null;) {
                if (doc.trim().length() == 0) {

                    System.out.println(numShortDocuments);
                    continue;
                }
                String[] words = doc.trim().split("\\s+");
                List<Integer> document = new ArrayList<Integer>();

                if(words.length==0)
                    System.out.println("here!");
                for (String word : words) {
                    if (word2IdVocabulary.containsKey(word)) {
                        document.add(word2IdVocabulary.get(word));
                    }
                    else {
                        indexWord += 1;
                        word2IdVocabulary.put(word, indexWord);
                        id2WordVocabulary.put(indexWord, word);
                        document.add(indexWord);
                    }
                }
                numShortDocuments++;
                numWordsInShortCorpus += document.size();
                ShortCorpus.add(document);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }


        br = null;
        try {
            br = new BufferedReader(new FileReader(pseudoCorpusPath));
            for (String doc; (doc = br.readLine()) != null;) {
                if (doc.trim().length() == 0) {
                    System.out.println(numPseudoDocuments);
                    continue;
                }
                String[] words = doc.trim().split("\\s+");
                List<Integer> document = new ArrayList<Integer>();

                if(words.length==0)
                    System.out.println("here!");
                for (String word : words) {
                    if (word2IdVocabulary.containsKey(word)) {
                        document.add(word2IdVocabulary.get(word));
                    }
                    else {
                        indexWord += 1;
                        word2IdVocabulary.put(word, indexWord);
                        id2WordVocabulary.put(indexWord, word);
                        document.add(indexWord);
                    }
                }

                numPseudoDocuments++;
                numWordsInPseudoCorpus += document.size();
                pseudoCorpus.add(document);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        V = word2IdVocabulary.size(); // vocabularySize = indexWord

        phi = new double[K][V];
        pz = new double[K];
        pdz = new double[numShortDocuments][K];

        nw = new int[V][K];
        nwsum = new int[K];
        npd = new double[numPseudoDocuments][K];
        npdsum = new int[numPseudoDocuments];
        nsd = new int[numShortDocuments][K];
        nsdsum = new int[numShortDocuments];
        sz = new int[numShortDocuments][];
        Pz = new int[numPseudoDocuments][];

        // initialise count variables.
        multiPros = new double[K];
        for (int i = 0; i < K; i++) {
            multiPros[i] = 1.0 / K;
        }

        System.out.println("Corpus size: " + numShortDocuments + " docs, "
                + numWordsInShortCorpus + " words");
        System.out.println("Pseudo Corpus size: " + numPseudoDocuments + " docs, "
                + numWordsInPseudoCorpus + " words");
        System.out.println("Vocabuary size: " + V);
        System.out.println("Number of topics: " + K);
        System.out.println("alpha: " + alpha);
        System.out.println("beta: " + beta);
        System.out.println("Number of sampling iterations: " + ITERATIONS);
        System.out.println("Number of top topical words: " + topWords);

        initialize();
    }

    /**
     * Randomly initialize topic assignments
     */
    public void initialize()
            throws IOException
    {
        System.out.println("Randomly initializing topic assignments ...");

        long startTime = System.currentTimeMillis();
        for (int i = 0; i < numPseudoDocuments; i++) {

            int docSize = pseudoCorpus.get(i).size();
            Pz[i] = new int[docSize];
            for (int j = 0; j < docSize; j++) {
                int topic = FuncUtils.nextDiscrete(multiPros); // Sample a topic
                // Increase counts
                npd[i][topic] += 1;
                nw[pseudoCorpus.get(i).get(j)][topic] += 1;
                npdsum[i] += 1;
                nwsum[topic] += 1;
                Pz[i][j] = topic;
            }

        }
        for (int i = 0; i < numShortDocuments; i++) {
            int docSize = ShortCorpus.get(i).size();
            sz[i] = new int[docSize];
            for (int j = 0; j < docSize; j++) {
                int topic = FuncUtils.nextDiscrete(npd[i]); // Sample a topic
                // Increase counts
                nsd[i][topic] += 1;
                nw[ShortCorpus.get(i).get(j)][topic] += 1;
                nsdsum[i] += 1;
                nwsum[topic] += 1;
                sz[i][j] = topic;
            }
        }
        initTime =System.currentTimeMillis()-startTime;
    }

    public void inference()
            throws IOException
    {

		writeDictionary();
        System.out.println("Running Gibbs sampling inference: ");

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[numShortDocuments][K];
            phisum = new double[K][V];
            numstats = 0;
        }
        System.out.println("Sampling " + ITERATIONS
                + " iterations with burn-in of " + BURN_IN + " (B/S="
                + THIN_INTERVAL + ").");

        long startTime = System.currentTimeMillis();
        int i;
        for (i = 0; i < ITERATIONS; i++) {
            // for all z_i
            for (int m = 0; m < Pz.length; m++) {
                for (int n = 0; n < Pz[m].length; n++) {

                    int topic = sampleFullConditionalForG(m, n);
                    Pz[m][n] = topic;
                }
            }
            for (int m = 0; m < sz.length; m++) {
                for (int n = 0; n < sz[m].length; n++) {

                    int topic = sampleFullConditionalForS(m, n);
                    sz[m][n] = topic;
                }
            }
            if ((i < BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("B");
                dispcol++;
            }
            // display progress
            if ((i > BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("S");
                dispcol++;
            }
            if (dispcol >= 100) {
                System.out.println();
                dispcol = 0;
            }
        }

        iterTime =System.currentTimeMillis()-startTime;
        expName = orgExpName;
        System.out.println("Writing output from the last sample ...");
        write();
        System.out.println("Sampling completed!");

    }


    /**
     * Sample a topic z_i from the full conditional distribution: p(z_i = j |
     * z_-i, w) = (n_-i,j(w_i) + beta)/(n_-i,j(.) + W * beta) * (n_-i,j(d_i) +
     * alpha)/(n_-i,.(d_i) + K * alpha)
     *
     * @param m
     *            document
     * @param n
     *            word
     */
    private int sampleFullConditionalForG(int m, int n) {

        // remove z_i from the count variables
        int topic = Pz[m][n];
        nw[pseudoCorpus.get(m).get(n)][topic]--;
        double[] original_npd = npd[m];
        npd[m][topic]--;
        nwsum[topic]--;
        npdsum[m]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        double ano;
        for (int k = 0; k < K; k++) {
            p[k] = ((nw[pseudoCorpus.get(m).get(n)][k] + beta) / (nwsum[k] + V * beta))
                    * (npd[m][k] + alpha) ;// (ngdsum[m] + nsdsum[m] + K * alpha)
            if(npd[m][k]!=0){
                ano = Math.pow(original_npd[k]*npdsum[m]/((npdsum[m]+1.)*npd[m][k]),nsd[m][k]);
                p[k]*=ano;
            }
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[K - 1];
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }
        // add newly estimated z_i to count variables
        nw[pseudoCorpus.get(m).get(n)][topic]++;
        npd[m][topic]++;
        nwsum[topic]++;
        npdsum[m]++;

        return topic;
    }

    private int sampleFullConditionalForS(int m, int n) {

        // remove z_i from the count variables
        int topic = sz[m][n];
        nw[ShortCorpus.get(m).get(n)][topic]--;
        nsd[m][topic]--;
        nwsum[topic]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (nw[ShortCorpus.get(m).get(n)][k] + beta) / (nwsum[k] + V * beta)
                    * npd[m][k];
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[K - 1];
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }
        // add newly estimated z_i to count variables
        nw[ShortCorpus.get(m).get(n)][topic]++;
        nsd[m][topic]++;
        nwsum[topic]++;

        return topic;
    }

    public void writeParameters()
            throws IOException
    {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".paras"));
        writer.write("-model" + "\t" + expName);
        writer.write("\n-ntopics" + "\t" + K);
        writer.write("\n-alpha" + "\t" + alpha);
        writer.write("\n-beta" + "\t" + beta);
        writer.write("\n-niters" + "\t" + ITERATIONS);
        writer.write("\n-twords" + "\t" + topWords);
        writer.write("\n-name" + "\t" + expName);

        writer.write("\n-initiation time" + "\t" + initTime);
        writer.write("\n-one iteration time" + "\t" + iterTime/ITERATIONS);
        writer.write("\n-total time" + "\t" + (initTime+iterTime));

        writer.close();
    }

    public void writeDictionary()
            throws IOException
    {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".vocabulary"));
        for (int id = 0; id < V; id++)
            writer.write(id2WordVocabulary.get(id) + " " + id + "\n");
        writer.close();
    }

    public void writeTopicAssignments()
            throws IOException
    {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".topicAssignments"));
        for (int dIndex = 0; dIndex < numShortDocuments; dIndex++) {
            int docSize = ShortCorpus.get(dIndex).size();
            for (int wIndex = 0; wIndex < docSize; wIndex++) {
                writer.write(sz[dIndex][wIndex] + " ");
            }
            writer.write("\n");
        }
        writer.close();
    }

    public void writeTopTopicalWords()
            throws IOException
    {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".topWords"));

        for (int tIndex = 0; tIndex < K; tIndex++) {
            //writer.write("Topic" + new Integer(tIndex) + ":");

            Map<Integer, Integer> wordCount = new TreeMap<Integer, Integer>();
            for (int wIndex = 0; wIndex < V; wIndex++) {
                wordCount.put(wIndex, nw[wIndex][tIndex]);
            }
            wordCount = FuncUtils.sortByValueDescending(wordCount);

            Set<Integer> mostLikelyWords = wordCount.keySet();
            int count = 0;
            for (Integer index : mostLikelyWords) {
                if (count < topWords) {
                    double pro = (nw[index][tIndex] + beta)
                            / (nwsum[tIndex] + beta*V);
                    pro = Math.round(pro * 1000000.0) / 1000000.0;
                    writer.write( id2WordVocabulary.get(index) + " ");
                    count += 1;
                }
                else {
                    writer.write("\n");
                    break;
                }
            }
        }
        writer.close();
    }

    /**
     * Retrieve estimated topic--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phi[][] is taken.
     *

     */
    public void writeTopicWordPros()
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".phi"));
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    double pro = phisum[k][w] / numstats;
                    writer.write(pro + " ");
                }
                writer.write("\n");
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    double pro = (nw[w][k] + beta) / (nwsum[k] + V * beta);
                    writer.write(pro + " ");
                }
                writer.write("\n");
            }
        }
        writer.close();

    }

    /**
     * Retrieve estimated document--topic associations. If sample lag > 0 then
     * the mean value of all sampled statistics for theta[][] is taken.
     *
     */
    public void compute_phi() {
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < V; w++) {
                phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
            }
        }
    }

    public void compute_pz() {
        double sum = 0.0;
        for (int i = 0; i < K; i++) {
            sum += nwsum[i];
        }
        for (int i = 0; i < K; i++) {
            pz[i] = (nwsum[i] + alpha) / (sum + alpha*K);
        }
    }

    public void compute_pzd() {
        double[][] pwz = new double[V][K]; // pwz[word][topic]
        for (int i = 0; i < V; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < K; j++) {
                pwz[i][j] = pz[j] * phi[j][i];
                row_sum += pwz[i][j];
            }
            for (int j = 0; j < K; j++) {
                pwz[i][j] = pwz[i][j] / row_sum;
            }

        }

        for (int i = 0; i < numShortDocuments; i++) {
            List<Integer> doc_word_id = ShortCorpus.get(i);
            double row_sum = 0.0;
            for (int j = 0; j < K; j++) {

                for (int wordID : doc_word_id) {
                    pdz[i][j] += pwz[wordID][j];
                }
                row_sum += pdz[i][j];

            }
            for (int j = 0; j < K; j++) {
                pdz[i][j] = pdz[i][j] / row_sum;
            }
        }
    }
    public void writeDocTopicPros()
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".theta"));
        for (int i = 0; i < numShortDocuments; i++) {
            for (int tIndex = 0; tIndex < K; tIndex++) {
                writer.write((pdz[i][tIndex]) + " ");
            }
            writer.write("\n");
        }
        writer.close();
    }

    public void write()
            throws IOException {
        compute_phi();
        compute_pz();
        compute_pzd();
        writeParameters();
        writeTopTopicalWords();
        writeDocTopicPros();
        writeTopicAssignments();
    }

    public static void main(String[] args)
            throws Exception
    {
        IETM ietm = new IETM("./dataset/Tweet/Tweet.txt",
                "./dataset/Tweet/Tweet_GPT.txt",
                "./results/Tweet/",50, 1.,
                0.01,1000, 20, "IETM_50_1.0_0.01_1");
        ietm.inference();
    }
}

