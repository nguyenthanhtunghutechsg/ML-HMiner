package MLHUI_Miner2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// ML-HUI-MINER ALGORITHM MULTI-CORE
// ---------------------------------
// This algorithm extends the FHM algorithm to mine multi-level high-utility itemsets
// 
// Coded by Trinh D.D. Nguyen
// Version 1.1 - Nov, 2020
//
public class AlgoMLHUIMinerMC {

	public long timerStart = 0; // time stamps for benchmarking purpose
	public long timerStop = 0;
	public long patternCount = 0; // multi-level HUIs found
	public long candidateCount = 0; // candidate high-utility itemsets counter
	public int transCount = 0; // number of transactions
	public double minUtil = 0.0; // minimum utility
	int itemCount[];
	Map<Integer, Integer> mapItemToLevel;
	double twus[]; // Map to remember the TWU of each item
	Map<Integer, List<Integer>> mapItemToAncestor;
	List<double[][]> EUCSPerLevel; // EUCS: key:item key:another item value:twu
	Dataset dataset;
	int BUFFERS_SIZE = 50000; // for storing the current itemset that is mined during the mining process
	int[] itemsetBuffer = null;
	ArrayList<List<UtilityList>> UtilityListPerLevel;
	ArrayList<int[]> oldNameToNewNamesPerLevel;
	ArrayList<int[]> newNamesToOldNamesPerLevel;
	Taxonomy taxonomy; // for describing the taxonomy of a dataset

	class Pair { // represents an item and its utility in a transaction
		int item = 0;
		double utility = 0.0;
	}

	// constructor
	public AlgoMLHUIMinerMC() {
	}

	// main ML-HUI-Miner MC algorithm
	public void runAlgorithm(String inputTransaction, String inputTaxonomy, String output, Double minUtility)
			throws IOException {

		// initializations
		MemoryLogger.getInstance().reset();
		itemsetBuffer = new int[BUFFERS_SIZE];
		EUCSPerLevel = new ArrayList<double[][]>();
		mapItemToLevel = new HashMap<Integer, Integer>();
		mapItemToAncestor = new HashMap<Integer, List<Integer>>();
		dataset = new Dataset(inputTransaction, Integer.MAX_VALUE);
		taxonomy = new Taxonomy(inputTaxonomy, dataset);
		minUtil = minUtility;
		UtilityListPerLevel = new ArrayList<>();
		timerStart = System.currentTimeMillis();

		// first dataset scan to calculate the TWU of each item.
		System.out.println("First dataset scan...");

		transCount = dataset.getTransactions().size();
		calculateTWU();
		int maxLevel = getMaxLevel(mapItemToLevel);
		ArrayList<ArrayList<Integer>> itemsToKeepPerLevel = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < maxLevel; i++) {
			itemsToKeepPerLevel.add(new ArrayList<Integer>());
			itemCount = new int[maxLevel];

		}
		for (int j = 1; j < twus.length; j++) {
			if (twus[j] >= minUtil) {
				itemsToKeepPerLevel.get(mapItemToLevel.get(j) - 1).add(j);
			}
		}
		insertionSort(itemsToKeepPerLevel, twus);
		oldNameToNewNamesPerLevel = new ArrayList<int[]>();
		newNamesToOldNamesPerLevel = new ArrayList<int[]>();

		for (int i = 0; i < maxLevel; i++) {
			ArrayList<Integer> itemsToKeep = itemsToKeepPerLevel.get(i);
			int ItemLevel = itemsToKeep.size();
			itemCount[i] = ItemLevel;
			double[][] EUCS = new double[ItemLevel + 1][ItemLevel + 1];
			EUCSPerLevel.add(EUCS);
			int[] oldNameToNewNames = new int[dataset.getMaxItem() + 1];
			// This structure will store the old name corresponding to each new name
			int[] newNamesToOldNames = new int[dataset.getMaxItem() + 1];
			int currentName = 1;

			// For each item in increasing order of TWU
			for (int j = 0; j < itemsToKeep.size(); j++) {
				// get the item old name
				int item = itemsToKeep.get(j);
				// give it the new name
				oldNameToNewNames[item] = currentName;
				// remember its old name
				newNamesToOldNames[currentName] = item;
				// replace its old name by the new name in the list of promising items
				itemsToKeep.set(j, currentName);
				// increment by one the current name so that
				currentName++;
			}
			oldNameToNewNamesPerLevel.add(oldNameToNewNames);
			newNamesToOldNamesPerLevel.add(newNamesToOldNames);
		}

		for (int i = 0; i < maxLevel; i++) {
			List<UtilityList> listOfUtilityLists = new ArrayList<UtilityList>();
			for (int j = 0; j < itemCount[i]; j++) {
				listOfUtilityLists.add(new UtilityList(j + 1));

			}
			UtilityListPerLevel.add(listOfUtilityLists);
		}
		for (int i = 0; i < dataset.getTransactions().size(); i++) {
			Transaction transaction = dataset.transactions.get(i);
			transaction.setLevelTransaction(maxLevel);
			transaction.removeUnpromisingItems(oldNameToNewNamesPerLevel, mapItemToAncestor, mapItemToLevel);
			// System.out.println(transaction.itemLevelToStringAll(0));
		}

		System.out.println("==== DATASET CHARACTERISTICS ====");
		System.out.println(" Transactions: " + transCount);
		System.out.println(" Levels      : " + getMaxLevel(mapItemToLevel));
		System.out.println(" |GI|        : " + taxonomy.parentCount());
		System.out.println("=================================");

		System.out.println("Construct.");
		int tid = 0;
		for (Transaction tran : dataset.transactions) {
			tid++;
			for (int i = 0; i < maxLevel; i++) {
				if (tran.listTransactionUtility.get(i) == 0) {
					continue;
				}
				double ru = 0;
				ArrayList<Integer> itemInTransactionInLevel = tran.ListItemPerLevel.get(i);
				ArrayList<Double> UtilityInTransactionInLevel = tran.ListUtilityPerLevel.get(i);
				for (int j = itemInTransactionInLevel.size() - 1; j >= 0; j--) {
					int item = itemInTransactionInLevel.get(j);
					double nU = UtilityInTransactionInLevel.get(j);
					Element element = new Element(tid, nU, ru);
					UtilityList utilityListOfItem = UtilityListPerLevel.get(i).get(item - 1);
					if (utilityListOfItem != null) {
						utilityListOfItem.addElement(element);
					}
					ru = ru + nU;
				}

			}
		}
		for (Transaction tran : dataset.transactions) {
			for (int l = 0; l < maxLevel; l++) {
				ArrayList<Integer> itemInTransactionInLevel = tran.ListItemPerLevel.get(l);
				int soluong = itemInTransactionInLevel.size();
				double tu = tran.listTransactionUtility.get(l);
				for (int i = 0; i < soluong - 1; i++)
					for (int j = i + 1; j < soluong; j++)
						EUCSPerLevel.get(l)[itemInTransactionInLevel.get(i)][itemInTransactionInLevel.get(j)] += tu;
			}

		}
		System.out.println("Mining.");
		for (int i = 1; i < maxLevel; i++) { // Mine the database recursively
			fhm(itemsetBuffer, 0, null, UtilityListPerLevel.get(i), minUtility, i);
		}
		MemoryLogger.getInstance().checkMemory(); // check the memory usage again and close the file.
		timerStop = System.currentTimeMillis(); // record end time

		System.out.println("Done.");
	}

	// Method to compare items by their TWU

	// FHM algorithm
	private void fhm(int[] prefix, int prefixLength, UtilityList pUL, List<UtilityList> ULs, Double minUtility,
			int level) throws IOException {

		// For each extension X of prefix P
		for (int i = 0; i < ULs.size(); i++) {
			UtilityList X = ULs.get(i);

			// If pX is a high utility itemset. we save the itemset: pX
			if (X.sumIutils >= minUtility) {
				writeOut(prefix, prefixLength, X.item, X.sumIutils); // output
			}

			// If the sum of the remaining utilities for pX
			// is higher than minUtility, we explore extensions of pX.
			// (this is the pruning condition)
			if (X.sumIutils + X.sumRutils >= minUtility) {
				// This list will contain the utility lists of pX extensions.
				List<UtilityList> exULs = new ArrayList<UtilityList>();
				// For each extension of p appearing
				// after X according to the ascending order
				for (int j = i + 1; j < ULs.size(); j++) {
					UtilityList Y = ULs.get(j);

					// ======================== NEW OPTIMIZATION USED IN FHM
					Double twuF = EUCSPerLevel.get(level)[X.item][Y.item];

					if (twuF < minUtility) {
						continue;
					}

					candidateCount++;
					// =========================== END OF NEW OPTIMIZATION

					// we construct the extension pXY
					// and add it to the list of extensions of pX
					UtilityList temp = construct(pUL, X, Y, minUtility);
					if (temp != null) {
						exULs.add(temp);
					}
				}

				// We create new prefix pX
				itemsetBuffer[prefixLength] = X.item;

				// We make a recursive call to discover all itemsets with the prefix pXY
				fhm(itemsetBuffer, prefixLength + 1, X, exULs, minUtility, level);
			}
		}
		MemoryLogger.getInstance().checkMemory();
	}

	// construct a utility list
	private UtilityList construct(UtilityList P, UtilityList px, UtilityList py, Double minUtility) {
		// create an empy utility list for pXY
		UtilityList pxyUL = new UtilityList(py.item);

		// == new optimization - LA-prune == /
		// Initialize the sum of total utility
		double totalUtility = px.sumIutils + px.sumRutils;
		// ================================================

		// for each element in the utility list of pX
		for (Element ex : px.elements) {
			// do a binary search to find element ey in py with tid = ex.tid
			Element ey = findElementWithTID(py, ex.tid);
			if (ey == null) {
				totalUtility -= (ex.iutils + ex.rutils);
				if (totalUtility < minUtility) {
					return null;
				}
				// =============================================== /
				continue;
			}
			// if the prefix p is null
			if (P == null) {
				// Create the new element
				Element eXY = new Element(ex.tid, ex.iutils + ey.iutils, ey.rutils);
				// add the new element to the utility list of pXY
				pxyUL.addElement(eXY);

			} else {
				// find the element in the utility list of p wih the same tid
				Element e = findElementWithTID(P, ex.tid);
				if (e != null) {
					// Create new element
					Element eXY = new Element(ex.tid, ex.iutils + ey.iutils - e.iutils, ey.rutils);
					// add the new element to the utility list of pXY
					pxyUL.addElement(eXY);
				}
			}
		}
		// return the utility list of pXY.
		return pxyUL;
	}

	// get maximum level from the taxonomy
	private static Integer getMaxLevel(Map<Integer, Integer> map) {
		if (map == null)
			return null;
		int length = map.size();
		Collection<Integer> c = map.values();
		Object[] obj = c.toArray();
		Arrays.sort(obj);
		return Integer.parseInt(obj[length - 1].toString());
	}

	// Do a binary search to find the element with a given tid in a utility list
	private Element findElementWithTID(UtilityList ulist, int tid) {
		List<Element> list = ulist.elements;

		int first = 0; // do a search to check if the subset appears in level k-1.
		int last = list.size() - 1;

		while (first <= last) { // the binary search
			int middle = (first + last) >>> 1;
			if (list.get(middle).tid < tid)
				first = middle + 1;
			else if (list.get(middle).tid > tid)
				last = middle - 1;
			else
				return list.get(middle);
		}
		return null;
	}

	// output the pattern
	// we bypass this step in multi-threaded environment, just increase the pattern
	// count
	private void writeOut(int[] prefix, int prefixLength, int item, double utility) throws IOException {
		patternCount++; // increase the number of high utility itemsets found
		StringBuilder buffer = new StringBuilder();
		// append the prefix
		for (int i = 0; i < prefixLength; i++) {
			buffer.append(newNamesToOldNamesPerLevel.get(1)[prefix[i]]);
			buffer.append(' ');
		}
		// append the last item
		buffer.append(newNamesToOldNamesPerLevel.get(1)[item]);
		// append the utility value
		buffer.append(" #UTIL: ");
		buffer.append(utility);
		// write to file
		System.out.println(buffer.toString());
	}

	public void calculateTWU() {

		// Initialize utility bins for all items
		twus = new double[dataset.getMaxItem() + 1];

		// Scan the database to fill the utility bins
		// For each transaction
		transCount = dataset.getTransactions().size();
		for (int tid = 0; tid < transCount; tid++) {

			Transaction transaction = dataset.getTransactions().get(tid);
			ArrayList<Integer> ancestantExist = new ArrayList<Integer>();

			for (int i = 0; i < transaction.getItems().length; i++) { // for each item, add the transaction utility to
																		// its TWU
				Integer item = transaction.getItems()[i];
				double transactionUtility = transaction.getTu();
				Double twu = twus[item]; // get the current TWU of that item

				// add the utility of the item in the current transaction to its twu
				twus[item] = twu + transactionUtility;

				ArrayList<Integer> ancestor = new ArrayList<Integer>();

				ancestor.add(item);
				if (mapItemToAncestor.get(item) == null) {
					Integer itemCopy = item;
//					for (int m = 0; m < taxonomy.size(); m++) {
//
//						Integer childItem = taxonomy.child(m);								
//						Integer parentItem = taxonomy.parent(m);
//						
//						if (childItem.intValue() == itemCopy.intValue()) {
//							ancestor.add(parentItem);
//							if (!ancestantExist.contains(parentItem)) {
//								ancestantExist.add(parentItem);
//								Double twuParent = twus[parentItem];
//								twuParent += transactionUtility;
//								twus[parentItem] = twuParent;
//							}
//							itemCopy = parentItem;
//						}
//					}
					while (itemCopy != null) {
						Integer childItem = itemCopy;
						Integer parentItem = taxonomy.MapdataParent.get(childItem);
						if (parentItem != null) {
							ancestor.add(parentItem);
							if (!ancestantExist.contains(parentItem)) {
								ancestantExist.add(parentItem);
								Double twuParent = twus[parentItem];
								twuParent = (twuParent == null) ? transactionUtility : transactionUtility + twuParent;
								twus[parentItem] = twuParent;
							}
						}
						itemCopy = parentItem;
					}
					int k = ancestor.size();
					for (int j = ancestor.size() - 1; j >= 0; j--, k--) {
						if (mapItemToLevel.get(ancestor.get(j)) == null) {
							mapItemToLevel.put(ancestor.get(j), k);
						} else {
							if (k < mapItemToLevel.get(ancestor.get(j))) {
								mapItemToLevel.put(ancestor.get(j), k);
							}

						}
					}
					for (int itemKey = 0; itemKey < ancestor.size(); itemKey++) {
						List<Integer> itemValue = new ArrayList<>();
						for (int listValue = itemKey; listValue < ancestor.size(); listValue++) {
							itemValue.add(ancestor.get(listValue));

						}
						mapItemToAncestor.put(ancestor.get(itemKey), itemValue);
					}

				} else {
					List<Integer> listAncestorOfItem = mapItemToAncestor.get(item);

					for (int k = 1; k < listAncestorOfItem.size(); k++) {
						if (!ancestantExist.contains(listAncestorOfItem.get(k))) {
							ancestantExist.add(listAncestorOfItem.get(k));
							Double twuParent = twus[listAncestorOfItem.get(k)];
							twuParent += transaction.getTu();
							twus[listAncestorOfItem.get(k)] = twuParent;
						}
					}
				}
			}
		}
	}

	// print statistics
	public void printStats() throws IOException {
		System.out.println("=============  ML-HUI-Miner MC STATISTICS =============");
		System.out.println(" Given minutil     : " + minUtil);
		System.out.println(" Approx runtime    : " + (timerStop - timerStart) + " ms");
		System.out.println(" Approx memory used: " + MemoryLogger.getInstance().getMaxMemory() + " MB");
		System.out.println(" MLHUIs found      : " + patternCount);
		System.out.println(" Candidate count   : " + candidateCount);
		System.out.println("=======================================================");
	}

	public static void insertionSort(ArrayList<ArrayList<Integer>> itemList, double[] ArrayTWU) {
		// the following lines are simply a modified an insertion sort
		for (List<Integer> items : itemList) {
			for (int j = 1; j < items.size(); j++) {
				Integer itemJ = items.get(j);
				int i = j - 1;
				Integer itemI = items.get(i);

				// we compare the twu of items i and j
				double comparison = ArrayTWU[itemI] - ArrayTWU[itemJ];
				// if the twu is equal, we use the lexicographical order to decide whether i is
				// greater
				// than j or not.
				if (comparison == 0) {
					comparison = itemI - itemJ;
				}

				while (comparison > 0) {
					items.set(i + 1, itemI);

					i--;
					if (i < 0) {
						break;
					}

					itemI = items.get(i);
					comparison = ArrayTWU[itemI] - ArrayTWU[itemJ];
					// if the twu is equal, we use the lexicographical order to decide whether i is
					// greater
					// than j or not.
					if (comparison == 0) {
						comparison = itemI - itemJ;
					}
				}
				items.set(i + 1, itemJ);
			}
		}
	}

}