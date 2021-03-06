package CLHMiner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class pCLH_Miner {
	double maxMemory = 0; 
	int minUtil;
	List<UtilityList> ListUls;
	int itemCount = 0;
	int giCount = 0;
	int taxDepth = 0;
	public long startTimestamp = 0;
	Map<Integer, Double> mapItemToTWU;
	/** the time at which the algorithm ended */
	public long endTimestamp = 0;
	TaxonomyTree taxonomy;
	List<Integer> listItemPromsing;
	static int countHUI;
	static int candidate;
	List<Integer> ListItemLevel1Global;

	class Pair {
		int item = 0;
		double utility = 0;
	}

	class TransactionPair {
		List<Pair> ListPair;
		double TU;

		public TransactionPair() {
			ListPair = new ArrayList<Pair>();
			TU = 0;
		}
	}

	public void runAlgorithm(int minUtil, String inputPath, String outputPath, String TaxonomyPath) throws IOException {
		this.minUtil = minUtil;
		candidate = 0;
		maxMemory =0;
		startTimestamp = System.currentTimeMillis();
		mapItemToTWU = new HashMap<Integer, Double>();
		taxonomy = new TaxonomyTree();
		taxonomy.ReadDataFromPath(TaxonomyPath);
		BufferedReader myInput = null;
		ListItemLevel1Global = new ArrayList<Integer>();
		List<TransactionPair> datasetAfterRemove = new ArrayList<>();
		countHUI = 0;
		Set<Integer> itemInDB = new HashSet<Integer>();
		String thisLine;
		try {
			// prepare the object for reading the file
			myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputPath))));
			// for each line (transaction) until the end of file
			while ((thisLine = myInput.readLine()) != null) {
				// if the line is a comment, is empty or is a
				// kind of metadata
				if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%'
						|| thisLine.charAt(0) == '@') {
					continue;
				}

				// split the transaction according to the : separator
				String split[] = thisLine.split(":");
				// the first part is the list of items
				String items[] = split[0].split(" ");
				// the second part is the transaction utility
				double transactionUtility = Double.parseDouble(split[1]);
				HashSet<Integer> setParent = new HashSet<Integer>();
				// for each item, we add the transaction utility to its TWU
				for (int i = 0; i < items.length; i++) {
					// convert item to integer
					Integer item = Integer.parseInt(items[i]);
					itemInDB.add(item);
					if (taxonomy.mapItemToTaxonomyNode.get(item) == null) {
						TaxonomyNode newNode = new TaxonomyNode(item);
						taxonomy.mapItemToTaxonomyNode.get(-1).addChildren(newNode);
						taxonomy.mapItemToTaxonomyNode.put(item, newNode);
					} else {
						TaxonomyNode parentNode = taxonomy.mapItemToTaxonomyNode.get(item).getParent();
						while (parentNode.getData() != -1) {
							setParent.add(parentNode.getData());
							parentNode = parentNode.getParent();
						}
					}

					// get the current TWU of that item
					Double twu = mapItemToTWU.get(item);
					// add the utility of the item in the current transaction to its twu
					twu = (twu == null) ? transactionUtility : twu + transactionUtility;
					mapItemToTWU.put(item, twu);
				}
				for (Integer parentItemInTransaction : setParent) {
					Double twu = mapItemToTWU.get(parentItemInTransaction);
					twu = (twu == null) ? transactionUtility : twu + transactionUtility;
					mapItemToTWU.put(parentItemInTransaction, twu);
				}
			}
		} catch (Exception e) {
			// catches exception if error while reading the input file
			e.printStackTrace();
		} finally {
			if (myInput != null) {
				myInput.close();
			}
		}
		listItemPromsing = new ArrayList<Integer>();

		// For each item
		for (Integer item : mapItemToTWU.keySet()) {
			// if the item is promising (TWU >= minutility)
			if (mapItemToTWU.get(item) >= minUtil) {
				listItemPromsing.add(item);

			}
		}

		Collections.sort(listItemPromsing, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				// compare the TWU of the items
				return compareItems(o1, o2);
			}
		});

		try {
			myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputPath))));

			while ((thisLine = myInput.readLine()) != null) {
				// if the line is a comment, is empty or is a
				// kind of metadata
				if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%'
						|| thisLine.charAt(0) == '@') {
					continue;
				}
				String split[] = thisLine.split(":");
				// get the list of items
				String items[] = split[0].split(" ");
				// get the list of utility values corresponding to each item
				// for that transaction
				String utilityValues[] = split[2].split(" ");

				// Copy the transaction into lists but
				// without items with TWU < minutility

				// long newTWU = 0; // NEW OPTIMIZATION
				double TU = Double.parseDouble(split[1]);
				// Create a list to store items
				List<Pair> revisedTransaction = new ArrayList<Pair>();

				for (int i = 0; i < items.length; i++) {
					Double Utiliy = Double.parseDouble(utilityValues[i]);
					int item = Integer.parseInt(items[i]);
					Pair pair = new Pair();
					pair.item = item;
					pair.utility = Utiliy;
					revisedTransaction.add(pair);
				}
				

				TransactionPair transactionPair = new TransactionPair();
				transactionPair.ListPair = revisedTransaction;
				transactionPair.TU = TU;
				datasetAfterRemove.add(transactionPair);
			}

			List<Integer> listItemLevel1 = new ArrayList<Integer>();
			for (Integer item : listItemPromsing) {
				if (taxonomy.getMapItemToTaxonomyNode().get(item).getLevel() == 1) {
					listItemLevel1.add(item);
					ListItemLevel1Global.add(item);
				}

				if (taxonomy.getMapItemToTaxonomyNode().get(item).getLevel() > 1) {
					break;
				}
			}

			itemCount = itemInDB.size();
			giCount = taxonomy.getGI() - 1;
			taxDepth = taxonomy.getMaxLevel();
			listItemLevel1.stream().parallel().forEach(s -> {
				InitUtitlitList(datasetAfterRemove, s.intValue());
			});

		} catch (Exception e) {
			// TODO: handle exception
			throw e;
		}
		checkMemory();
		endTimestamp = System.currentTimeMillis();
		myInput.close();

	}

	private void InitUtitlitList(List<TransactionPair> Dataset, int itemID) {
		List<TransactionPair> NewDatasetAfterRemove = new ArrayList<pCLH_Miner.TransactionPair>();
		Map<Integer, UtilityList> mapItemToUtilityList = new HashMap<Integer, UtilityList>();
		for (int tid = 0; tid < Dataset.size(); tid++) {
			double remainingUtility = 0;
			List<Pair> revisedTransaction = Dataset.get(tid).ListPair;
			List<Pair> NewRevisedTransaction = new ArrayList<pCLH_Miner.Pair>();
			HashMap<Integer, Double> mapParentToUtility = new HashMap<Integer, Double>();
			for (int i = 0; i < revisedTransaction.size(); i++) {

				int item = revisedTransaction.get(i).item;
				double Utiliy = revisedTransaction.get(i).utility;
				TaxonomyNode nodeParent = taxonomy.mapItemToTaxonomyNode.get(item).getParent();
				while (nodeParent.getData() != -1) {
					Double utilityOfParent = mapParentToUtility.get(nodeParent.getData());
					if (utilityOfParent != null) {
						mapParentToUtility.put(nodeParent.getData(), utilityOfParent + Utiliy);
					} else {
						mapParentToUtility.put(nodeParent.getData(), Utiliy);
					}
					nodeParent = nodeParent.getParent();
				}
				Pair pair = new Pair();
				pair.item = item;
				pair.utility = Utiliy;
				if (mapItemToTWU.get(item) >= minUtil) {
					NewRevisedTransaction.add(pair);
					remainingUtility += pair.utility;
				}
			}
			Collections.sort(NewRevisedTransaction, new Comparator<Pair>() {
				public int compare(Pair o1, Pair o2) {
					return compareItems(o1.item, o2.item);
				}
			});
			double CountUtility = remainingUtility;
			for (int i = 0; i < NewRevisedTransaction.size(); i++) {
				Pair pair = NewRevisedTransaction.get(i);
				remainingUtility = remainingUtility - pair.utility;
				Element element = new Element(tid, pair.utility, remainingUtility, Dataset.get(tid).TU);
				UtilityList utilityListOfItem = mapItemToUtilityList.get(pair.item);
				if (utilityListOfItem == null) {
					UtilityList newUl = new UtilityList(pair.item);
					newUl.addElement(element);
					mapItemToUtilityList.put(pair.item, newUl);

				} else {
					utilityListOfItem.addElement(element);
				}

			}
			for (Integer itemParent : mapParentToUtility.keySet()) {
				double CountUtilityOfEachItem = CountUtility;
				for (int i = 0; i < NewRevisedTransaction.size(); i++) {
					Pair CurrentItem = NewRevisedTransaction.get(i);
					if (CheckParent(itemParent, CurrentItem.item)) {
						CountUtilityOfEachItem -= CurrentItem.utility;
					} else {
						if (compareItems(itemParent, CurrentItem.item) > 0) {
							CountUtilityOfEachItem -= CurrentItem.utility;
						}
					}
				}
				UtilityList utilityListOfItem = mapItemToUtilityList.get(itemParent);
				Element element = new Element(tid, mapParentToUtility.get(itemParent), CountUtilityOfEachItem,
						Dataset.get(tid).TU);
				if (utilityListOfItem != null) {
					utilityListOfItem.addElement(element);
				} else {
					UtilityList newUl = new UtilityList(itemParent);
					newUl.addElement(element);
					mapItemToUtilityList.put(itemParent, newUl);
				}

			}
			TransactionPair transactionPair = new TransactionPair();
			transactionPair.ListPair = NewRevisedTransaction;
			transactionPair.TU = Dataset.get(tid).TU;
			NewDatasetAfterRemove.add(transactionPair);

		}

		List<UtilityList> listUL = new ArrayList<UtilityList>();

		for (Integer item : ListItemLevel1Global) {
			listUL.add(mapItemToUtilityList.get(item));
		}
		
		
		SearchTree(new int[0], null, listUL, mapItemToUtilityList, NewDatasetAfterRemove,itemID);
	}

	private void SearchTree(int[] prefix, UtilityList pUL, List<UtilityList> ULs,
			Map<Integer, UtilityList> mapItemToUtilityList, List<TransactionPair> NewDatasetAfterRemove,int itemID) {
		for (int i = 0; i < ULs.size(); i++) {
			UtilityList X = ULs.get(i);
			if (X.sumIutils + X.sumRutils > minUtil) {
				TaxonomyNode taxonomyNodeX = taxonomy.getMapItemToTaxonomyNode().get(X.item);
				List<TaxonomyNode> childOfX = taxonomyNodeX.getChildren();
				for (TaxonomyNode taxonomyNode : childOfX) {
					int Child = taxonomyNode.getData();
					UtilityList ULofChild = mapItemToUtilityList.get(Child);
					if (ULofChild != null) {
						UtilityList exULBuild = constructTax(pUL, ULofChild, NewDatasetAfterRemove);
						X.AddChild(exULBuild);
					}
				}
				for (UtilityList childULs : X.getChild()) {
					if (childULs.GWU > minUtil) {
						ULs.add(childULs);
					}
				}
			}
			if (prefix.length==0) {
				if (X.item==itemID||CheckParent(itemID, X.item)) {
				//if(true) {
					candidate++;
					if (X.sumIutils > minUtil) {
						countHUI++;
						/*String result = "";
						for (int j = 0; j < prefix.length; j++) {
							result = result + " " + prefix[j];
						}
						result += " " + X.item + " ";
						// result += " #UTIL: " + X.sumIutils;
						System.out.println(result);*/

					}
					List<UtilityList> exULs = new ArrayList<UtilityList>();
					for (int j = i + 1; j < ULs.size(); j++) {

						UtilityList Y = ULs.get(j);

						if (!CheckParent(Y.item, X.item)) {
							UtilityList exULBuild = construct(pUL, X, Y, NewDatasetAfterRemove);
							if (exULBuild.GWU > minUtil) {
								exULs.add(exULBuild);
							}
						}

					}
					int[] newPrefix = new int[prefix.length + 1];
					System.arraycopy(prefix, 0, newPrefix, 0, prefix.length);
					newPrefix[prefix.length] = X.item;
					SearchTree(newPrefix, X, exULs, mapItemToUtilityList, NewDatasetAfterRemove,itemID);
				}
				
			}		
			else {
				candidate++;
				if (X.sumIutils > minUtil) {
					countHUI++;
					/*String result = "";
					for (int j = 0; j < prefix.length; j++) {
						result = result + " " + prefix[j];
					}
					result += " " + X.item + " ";
					// result += " #UTIL: " + X.sumIutils;
					System.out.println(result);*/

				}
				List<UtilityList> exULs = new ArrayList<UtilityList>();
				for (int j = i + 1; j < ULs.size(); j++) {

					UtilityList Y = ULs.get(j);

					if (!CheckParent(Y.item, X.item)) {
						UtilityList exULBuild = construct(pUL, X, Y, NewDatasetAfterRemove);
						if (exULBuild.GWU > minUtil) {
							exULs.add(exULBuild);
						}
					}

				}
				int[] newPrefix = new int[prefix.length + 1];
				System.arraycopy(prefix, 0, newPrefix, 0, prefix.length);
				newPrefix[prefix.length] = X.item;
				SearchTree(newPrefix, X, exULs, mapItemToUtilityList, NewDatasetAfterRemove,itemID);
			}
			
			
		}

	}

	private UtilityList constructTax(UtilityList P, UtilityList Child, List<TransactionPair> NewDatasetAfterRemove) {

		if (P == null) {
			return Child;
		} else {
			UtilityList newULs = new UtilityList(Child.item);
			for (Element PElment : P.getElement()) {

				Element UnionChild = findElementWithTID(Child, PElment.tid);
				if (UnionChild != null) {
					List<Pair> trans = NewDatasetAfterRemove.get(UnionChild.tid).ListPair;
					double remainUtility = 0;
					for (int i = 0; i < trans.size(); i++) {
						Integer currentItem = trans.get(i).item;
						if (compareItems(currentItem, Child.item) > 0 && (!CheckParent(Child.item, currentItem))
								&& (!CheckParent(Child.item, currentItem))) {
							remainUtility += trans.get(i).utility;
						}
					}

					// Create new element
					Element newElment = new Element(UnionChild.tid, PElment.iutils + UnionChild.iutils, remainUtility,
							UnionChild.TU);
					// add the new element to the utility list of pXY
					newULs.addElement(newElment);
				}
			}
			// return the utility list of pXY.
			return newULs;
		}
	}

	private UtilityList construct(UtilityList P, UtilityList px, UtilityList py,
			List<TransactionPair> NewDatasetAfterRemove) {
		UtilityList pxyUL = new UtilityList(py.item);

		// for each element in the utility list of pX
		for (Element ex : px.elements) {
			// do a binary search to find element ey in py with tid = ex.tid
			Element ey = findElementWithTID(py, ex.tid);
			if (ey == null) {
				continue;
			}
			// if the prefix p is null
			if (P == null) {
				// Create the new element
				List<Pair> trans = NewDatasetAfterRemove.get(ex.tid).ListPair;
				double remainUtility = 0;
				for (int i = 0; i < trans.size(); i++) {
					Integer currentItem = trans.get(i).item;
					if (compareItems(currentItem, py.item) > 0 && (!CheckParent(px.item, currentItem))
							&& (!CheckParent(py.item, currentItem))) {
						remainUtility += trans.get(i).utility;
					}
				}
				Element eXY = new Element(ex.tid, ex.iutils + ey.iutils, remainUtility, ey.TU);
				// add the new element to the utility list of pXY
				pxyUL.addElement(eXY);

			} else {
				// find the element in the utility list of p wih the same tid
				Element e = findElementWithTID(P, ex.tid);
				if (e != null) {
					List<Pair> trans = NewDatasetAfterRemove.get(e.tid).ListPair;
					double remainUtility = 0;
					for (int i = 0; i < trans.size(); i++) {
						Integer currentItem = trans.get(i).item;
						if (compareItems(currentItem, py.item) > 0 && (!CheckParent(px.item, currentItem))
								&& (!CheckParent(py.item, currentItem))) {
							remainUtility += trans.get(i).utility;
						}
					}
					Element eXY = new Element(ex.tid, ex.iutils + ey.iutils - e.iutils, remainUtility, ey.TU);
					// add the new element to the utility list of pXY
					pxyUL.addElement(eXY);
				}
			}
		}
		// return the utility list of pXY.
		return pxyUL;
	}

	private Element findElementWithTID(UtilityList ulist, int tid) {
		List<Element> list = ulist.elements;

		// perform a binary search to check if the subset appears in level k-1.
		int first = 0;
		int last = list.size() - 1;

		// the binary search
		while (first <= last) {
			int middle = (first + last) >>> 1; // divide by 2

			if (list.get(middle).tid < tid) {
				first = middle + 1; // the itemset compared is larger than the subset according to the lexical order
			} else if (list.get(middle).tid > tid) {
				last = middle - 1; // the itemset compared is smaller than the subset is smaller according to the
									// lexical order
			} else {
				return list.get(middle);
			}
		}
		return null;
	}

	private int compareItems(int item1, int item2) {
		int levelOfItem1 = taxonomy.getMapItemToTaxonomyNode().get(item1).getLevel();
		int levelOfItem2 = taxonomy.getMapItemToTaxonomyNode().get(item2).getLevel();
		if (levelOfItem1 == levelOfItem2) {
			int compare = (int) (mapItemToTWU.get(item1) - mapItemToTWU.get(item2));
			// if the same, use the lexical order otherwise use the TWU
			return (compare == 0) ? item1 - item2 : compare;
		} else {
			return levelOfItem1 - levelOfItem2;
		}
	}

	private boolean CheckParent(int item1, int item2) {
		TaxonomyNode nodeItem1 = taxonomy.getMapItemToTaxonomyNode().get(item1);
		TaxonomyNode nodeItem2 = taxonomy.getMapItemToTaxonomyNode().get(item2);
		int levelOfItem1 = nodeItem1.getLevel();
		int levelOfItem2 = nodeItem2.getLevel();
		if (levelOfItem1 == levelOfItem2) {
			return false;
		} else {
			if (levelOfItem1 > levelOfItem2) {
				TaxonomyNode parentItem1 = nodeItem1.getParent();
				while (parentItem1.getData() != -1) {
					if (parentItem1.getData() == nodeItem2.getData()) {
						return true;
					}
					parentItem1 = parentItem1.getParent();
				}
				return false;
			} else {
				TaxonomyNode parentItem2 = nodeItem2.getParent();
				while (parentItem2.getData() != -1) {
					if (parentItem2.getData() == nodeItem1.getData()) {
						return true;
					}
					parentItem2 = parentItem2.getParent();
				}
				return false;
			}
		}
	}

	public void printStats() throws IOException {
		System.out.println("=============  CLH-Miner =============");
		System.out.println(" |I|              : " + itemCount);
		System.out.println(" |GI|             : " + giCount);
		System.out.println(" Depth            : " + taxDepth);
		System.out.println(" Total time ~     : " + (endTimestamp - startTimestamp) + " ms");
		 System.out.println(" Memory ~ " + maxMemory +
		 " MB");
		System.out.println(" CLHUIs/Candidates: " + countHUI + "/" + candidate);
		System.out.println("======================================");
	}
	private void checkMemory() {
		// get the current memory usage
		double currentMemory = (Runtime.getRuntime().totalMemory() -  Runtime.getRuntime().freeMemory())
				/ 1024d / 1024d;
		// if higher than the maximum until now
		if (currentMemory > maxMemory) {
			// replace the maximum with the current memory usage
			maxMemory = currentMemory;
		}
	}
}