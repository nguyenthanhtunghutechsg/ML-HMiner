package CLHMiner;

import java.io.IOException;

public class Main {

	public static void main(String[] args) throws IOException {

//		String TaxonomyPath = "connectTaxonomy.txt";
//		String inputPath = "connect.txt";
		String TaxonomyPath = "fruithutTaxonomy.txt";
		String inputPath = "fruithut.txt";
		double minutil = 400000;
		for (int i = 0; i < 4; i++) {
			System.gc();
			CLH_Miner cl = new CLH_Miner();
			cl.runAlgorithm((int) minutil, inputPath, "", TaxonomyPath);

			System.out.println("minU: " + minutil);

			cl.printStats();
			minutil -= 100000;
		}

//2088282/2150177
	}
}
