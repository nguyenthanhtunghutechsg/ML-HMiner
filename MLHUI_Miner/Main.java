package MLHUI_Miner;

import java.io.IOException;

import MLHUI_Miner.AlgoMLHUIMinerMC;

// TEST DRIVE FOR ML-HUI-MINER MC
// Coded by Trinh D.D. Nguyen
public class Main {

	public static void main(String[] args) throws IOException {

		String dataset = "T10I4D100K";
		String trans = dataset + ".txt";
		String taxonomy = dataset + "Taxonomy.txt";
		double minutil = 5000;
		for (int i = 0; i < 10; i++) {
			AlgoMLHUIMinerMC mlhuiminer = new AlgoMLHUIMinerMC();
			mlhuiminer.runAlgorithm(trans, taxonomy, null, minutil);
			mlhuiminer.printStats();
			//minutil-=10000;
		}
	}

}
