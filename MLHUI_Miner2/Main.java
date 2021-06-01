package MLHUI_Miner2;

import java.io.IOException;



// TEST DRIVE FOR ML-HUI-MINER MC
// Coded by Trinh D.D. Nguyen
public class Main {

	public static void main(String[] args) throws IOException {

		String dataset = "fruithut";
		String trans = dataset + ".txt";
		String taxonomy = dataset + "Taxonomy.txt";
		double minutil =150000000;
		for (int i = 0; i < 1; i++) {
			System.gc();
			AlgoMLHUIMinerMC mlhuiminer = new AlgoMLHUIMinerMC();
			mlhuiminer.runAlgorithm(trans, taxonomy, null, minutil);
			mlhuiminer.printStats();
			minutil-=100000;
		}
	}

}
