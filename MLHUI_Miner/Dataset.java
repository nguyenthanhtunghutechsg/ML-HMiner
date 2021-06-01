package MLHUI_Miner;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

// ML-HUI-MINER ALGORITHM MULTI-CORE
// ---------------------------------
// This algorithm extends the FHM algorithm to mine multi-level high-utility itemsets
//
// Coded by Trinh D.D. Nguyen
// Version 1.1 - Nov, 2020
//

//This class represents a quantitative database
public class Dataset {
	
	List<Transaction> transactions;						// the list of transactions in this dataset
	private int maxItem = 0;							// the largest item name

	// main constructor
    public Dataset(String datasetPath, int maximumTransactionCount) throws IOException {

        transactions = new ArrayList<Transaction>();	// Initialize a list to store transactions in memory
  
        BufferedReader br = new BufferedReader(new FileReader(datasetPath));
        String line;
        int i = 0;

        while((line = br.readLine()) != null) { 

			if (line.isEmpty() == true || 				// bypass comments and empty lines
				line.charAt(0) == '#' ||  
				line.charAt(0) == '%' || 
				line.charAt(0) == '@') {
				continue;
			}
			i++;
			
			transactions.add(createTransaction(line));	// read the transaction
	
        	if(i==maximumTransactionCount) {			// prevent exceeding the number of transactions
        		break;
        	}			
        }
        br.close();
    }

    // create a transaction object from a string read from the input file
    private Transaction createTransaction(String line) {
    	
    	String[] split = line.split(":");								// split the line into tokens using ":"    	
    	double transactionUtility = Double.parseDouble(split[1]);		// Get the transaction utility
        String[] itemsString = split[0].split(" ");						// Get the list of items
        String[] itemsUtilitiesString = split[2].split(" ");			// Get the list of item utilities
        int[] items = new  int[itemsString.length];						// store the items and their utilities
        double[] utilities = new  double[itemsString.length];
        
        for (int i = 0; i < items.length; i++) {						// for each item        	
        	items[i] = Integer.parseInt(itemsString[i]);				// store that item        	
        	utilities[i] = Double.parseDouble(itemsUtilitiesString[i]);	// and its utility in that transaction
            
            // determine the largest item name 
            if(items[i] > maxItem) {
                maxItem = items[i];
            }
        }		
		return new Transaction(items, utilities, transactionUtility);	// create the transaction 
    }

    // returns the list of all transactions
    public List<Transaction> getTransactions() {
        return transactions;
    }

    // returns the largest item within the database
    public int getMaxItem() {
        return maxItem;
    }

}
