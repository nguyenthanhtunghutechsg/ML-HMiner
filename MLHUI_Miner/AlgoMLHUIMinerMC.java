package MLHUI_Miner;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import MLHUI_Miner.Dataset;
import MLHUI_Miner.MemoryLogger;
import MLHUI_Miner.Transaction;

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

	Map<Integer, Integer> mapItemToLevel;
	Map<Integer, Double> mapItemToTWU; // Map to remember the TWU of each item
	Map<Integer, List<Integer>> mapItemToAncestor;
	Map<Integer, Map<Integer, Double>> mapFMAP; // EUCS: key:item key:another item value:twu

	int BUFFERS_SIZE = 50000; // for storing the current itemset that is mined during the mining process
	int[] itemsetBuffer = null;

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
		mapFMAP = new HashMap<Integer, Map<Integer, Double>>();
		mapItemToTWU = new HashMap<Integer, Double>();
		mapItemToLevel = new HashMap<Integer, Integer>();
		mapItemToAncestor = new HashMap<Integer, List<Integer>>();
		taxonomy = new Taxonomy(inputTaxonomy);
		minUtil = minUtility;

		timerStart = System.currentTimeMillis();

		// first dataset scan to calculate the TWU of each item.
		System.out.println("First dataset scan...");

		Dataset dataset = new Dataset(inputTransaction, Integer.MAX_VALUE);
		transCount = dataset.getTransactions().size();

		for (int tid = 0; tid < transCount; tid++) {
			Transaction transaction = dataset.getTransactions().get(tid);
			ArrayList<Integer> ancestantExist = new ArrayList<Integer>();

			for (int i = 0; i < transaction.getItems().length; i++) { // for each item, add the transaction utility to
																		// its TWU
				Integer item = transaction.getItems()[i];
				double transactionUtility = transaction.getUtility();
				Double twu = mapItemToTWU.get(item); // get the current TWU of that item

				// add the utility of the item in the current transaction to its twu
				twu = (twu == null) ? transactionUtility : twu + transactionUtility;

				ArrayList<Integer> ancestor = new ArrayList<Integer>();

				ancestor.add(item);
				mapItemToTWU.put(item, twu);
				if (mapItemToAncestor.get(item) == null) {
					Integer itemCopy = item;
					while (itemCopy != null) {
						Integer childItem = itemCopy;
						Integer parentItem = taxonomy.MapdataParent.get(childItem);
						if (parentItem != null) {
							ancestor.add(parentItem);
							if (!ancestantExist.contains(parentItem)) {
								ancestantExist.add(parentItem);
								Double twuParent = mapItemToTWU.get(parentItem);
								twuParent = (twuParent == null) ? transactionUtility : transactionUtility + twuParent;
								mapItemToTWU.put(parentItem, twuParent);
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
							Double twuParent = mapItemToTWU.get(listAncestorOfItem.get(k));
							twuParent = (twuParent == null) ? transaction.getUtility()
									: twuParent + transaction.getUtility();
							mapItemToTWU.put(listAncestorOfItem.get(k), twuParent);
						}
					}
				}
			}
		}

//		for (Integer ints : mapItemToAncestor.keySet()) {
//			String s = "";
//			s = s + ints + " ";
//			for (Integer list : mapItemToAncestor.get(ints)) {
//				s = s + list + " ";
//			}
//			System.out.println(s);
//		}
		List<List<UtilityList>> listOfUtilityLists = new ArrayList<>(); // A LIST TO STORE THE UTILITY LIST OF ITEMS
																		// WITH TWU >= MIN_UTILITY.
		// CREATE A MAP TO STORE THE UTILITY LIST FOR EACH ITEM.
		// Key : item Value : utility list associated to that item
		Map<Integer, UtilityList> mapItemToUtilityList = new HashMap<Integer, UtilityList>();
		for (Integer item : mapItemToTWU.keySet()) { // For each item

			if (mapItemToTWU.get(item) >= minUtility) { // if the item is promising (TWU >= minutility)
				UtilityList uList = new UtilityList(item); // create an empty Utility List that we will fill later.
				mapItemToUtilityList.put(item, uList); // add the item to the list of high TWU items
			} else {
				List<Integer> listAncestorOfItem = mapItemToAncestor.get(item);
				for (int k = 0; k < listAncestorOfItem.size(); k++) {
					if (mapItemToTWU.get(listAncestorOfItem.get(k)) >= minUtility) {
						List<Integer> itemList = new ArrayList<Integer>();
						itemList.add(item);
						UtilityList tuList = new UtilityList(item);
						mapItemToUtilityList.put(item, tuList);
						break;
					}
				}

			}
		}

		List<List<List<Pair>>> revisedTransaction = new ArrayList<>();

		for (int i = 0; i < getMaxLevel(mapItemToLevel); i++) {

			List<List<Pair>> revisedTransactionTemp = new ArrayList<>();

			for (int j = 0; j < transCount; j++) {
				List<Pair> rrTemp = new ArrayList<Pair>();
				revisedTransactionTemp.add(rrTemp);
			}

			revisedTransaction.add(revisedTransactionTemp);
		}

		System.out.println("==== DATASET CHARACTERISTICS ====");
		System.out.println(" Transactions: " + transCount);
		System.out.println(" Levels      : " + getMaxLevel(mapItemToLevel));
		System.out.println(" |GI|        : " + taxonomy.parentCount());
		System.out.println("=================================");

		System.out.println("Second dataset scan...");
		int maxLevel = getMaxLevel(mapItemToLevel);

		for (int tid = 0; tid < transCount; tid++) {
			
			Transaction transaction = dataset.getTransactions().get(tid);
			int[] items = transaction.getItems();
			double[] remainingUtility = new double[maxLevel];
			double[] newTWU = new double[maxLevel];
			Map<Integer, Double> mapItemToUtility = new HashMap<Integer, Double>();
			double[] utilities = transaction.getUtilities();
			for (int i = 0; i < items.length; i++) {

				Integer item = transaction.getItems()[i];
				mapItemToUtility.put(item, utilities[i]);
				List<Integer> listParent = mapItemToAncestor.get(item);

				for (int k = 1; k < listParent.size(); k++) {
					int parentItem = listParent.get(k);
					Double UtilityOfParent = mapItemToUtility.get(parentItem);
					if (UtilityOfParent == null) {
						UtilityOfParent = utilities[i];
					} else {
						UtilityOfParent += utilities[i];
					}
					mapItemToUtility.put(parentItem, UtilityOfParent);
				}
				for (int j : mapItemToUtility.keySet()) {
					int level = mapItemToLevel.get(j);
					if (mapItemToTWU.get(j) > minUtil) {
						Pair pair = new Pair();
						pair.item = j;
						pair.utility = mapItemToUtility.get(j);
						revisedTransaction.get(level - 1).get(tid).add(pair);
						remainingUtility[level - 1] += pair.utility;
						newTWU[level - 1] += pair.utility;
					} else {

					}
				}

			}

			// sort the transaction
			for (int i = 0; i < /* getMaxLevel(mapItemToLevel) */maxLevel; i++) {
				Collections.sort(revisedTransaction.get(i).get(tid), new Comparator<Pair>() {
					public int compare(Pair o1, Pair o2) {
						return compareItems(o1.item, o2.item);
					}
				});
			}
			long Timesss = 0;
			for (int levels = /* getMaxLevel(mapItemToLevel) */maxLevel - 1; levels >= 0; levels--) {
				for (int i = 0; i < revisedTransaction.get(levels).get(tid).size(); i++) {
					Pair pair = revisedTransaction.get(levels).get(tid).get(i);

					// subtract the utility of this item from the remaining utility
					remainingUtility[levels] = remainingUtility[levels] - pair.utility;

					// get the utility list of this item
					UtilityList utilityListOfItem = mapItemToUtilityList.get(pair.item);

					// Add a new Element to the utility list of this item corresponding to this
					// transaction
					Element element = new Element(tid, pair.utility, remainingUtility[levels]);

					if (utilityListOfItem != null) {
						utilityListOfItem.addElement(element);
					}
					
					// BEGIN NEW OPTIMIZATION for FHM
//					Map<Integer, Double> mapFMAPItem = mapFMAP.get(pair.item);
//					if (mapFMAPItem == null) {
//						mapFMAPItem = new HashMap<Integer, Double>();
//						mapFMAP.put(pair.item, mapFMAPItem);
//					}
//					long times = System.currentTimeMillis();
//					for (int j = i + 1; j < revisedTransaction.get(levels).get(tid).size(); j++) {
//						Pair pairAfter = revisedTransaction.get(levels).get(tid).get(j);
//
//						Double twuSum = mapFMAPItem.get(pairAfter.item);
//						if (twuSum == null) {
//							mapFMAPItem.put(pairAfter.item, newTWU[levels]);
//						} else {
//							mapFMAPItem.put(pairAfter.item, twuSum + newTWU[levels]);
//						}
//					}
//					long timess = System.currentTimeMillis();
//					Timesss+=(timess-times);
					// END OPTIMIZATION of FHM
				}
			}
			//System.out.println(tid + "	" + Timesss);
		}

		System.out.println("Constructing utility lists...");
		for (int i = 0; i < maxLevel; i++) {

			List<UtilityList> UtilityListOfILevel = new ArrayList<>();
			for (Integer item : mapItemToTWU.keySet()) {
				if (mapItemToTWU.get(item) >= minUtility) { // if the item is promising (TWU >= minutility)
					if (mapItemToLevel.get(item) == i + 1) {
						UtilityList uList = mapItemToUtilityList.get(item); // create an empty Utility List that we will
																			// fill later.
						// add the item to the list of high TWU items
						UtilityListOfILevel.add(uList);
					}
				}
			}

			listOfUtilityLists.add(UtilityListOfILevel);

			// SORT THE LIST OF HIGH TWU ITEMS IN ASCENDING ORDER
			Collections.sort(listOfUtilityLists.get(i), new Comparator<UtilityList>() {
				public int compare(UtilityList o1, UtilityList o2) {
					return compareItems(o1.item, o2.item); // compare the TWU of the items
				}
			});
		}

		System.out.println("Mining...");
		MemoryLogger.getInstance().checkMemory(); // check the memory usage

		/* ==== SEQUENTIAL CODE HERE */
		for (int i = 0; i < maxLevel; i++) { // Mine the database recursively
			fhm(itemsetBuffer, 0, null, listOfUtilityLists.get(i), minUtility);
		}

		MemoryLogger.getInstance().checkMemory(); // check the memory usage again and close the file.
		timerStop = System.currentTimeMillis(); // record end time

		System.out.println("Done.");
	}

	// Method to compare items by their TWU
	private int compareItems(int item1, int item2) {
		int compare = (int) (mapItemToTWU.get(item1) - mapItemToTWU.get(item2));
		// if the same, use the lexical order otherwise use the TWU
		return (compare == 0) ? item1 - item2 : compare;
	}

	// FHM algorithm
	private void fhm(int[] prefix, int prefixLength, UtilityList pUL, List<UtilityList> ULs, Double minUtility)
			throws IOException {

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
//					Map<Integer, Double> mapTWUF = mapFMAP.get(X.item);
//					if (mapTWUF != null) {
//						Double twuF = mapTWUF.get(Y.item);
//						if (twuF == null || twuF < minUtility) {
//							continue;
//						}
//					}
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
				fhm(itemsetBuffer, prefixLength + 1, X, exULs, minUtility);
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
			buffer.append(prefix[i]);
			buffer.append(' ');
		}
		// append the last item
		buffer.append(item);
		// append the utility value
		buffer.append(" #UTIL: ");
		buffer.append(utility);
		// write to file
		// System.out.println(buffer.toString());
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

}