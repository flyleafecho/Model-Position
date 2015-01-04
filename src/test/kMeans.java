package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * k-means clustering algorithm
 * 
 * @author Antonio Severien
 * 
 */

class Point {
	List<Double> values;
	int cluster;
	int id;
	boolean isVirtual;
	Point(){
		values = new ArrayList<Double>();
	}
}
public class kMeans {

	File dataset;
	int kClusters;
	int dimension;
	int dataSetSize;
	List<Point> points;
	List<Point> clustersCentroids;
	List<List<Integer>> clustersMapping;

	/*public kMeans() {
		kClusters = 3;
		dimension = 2;
		points = new ArrayList<List<Double>>(dataSetSize);
		clustersCentroids = new ArrayList<List<Double>>(kClusters);
		clustersMapping = new ArrayList<List<Integer>>(kClusters);
		for (int i = 0; i < kClusters; i++) {
			clustersMapping.add(i, new ArrayList<Integer>());
		}
	}*/
	
	public kMeans(int clusters, int dimension, int dataSetSize) {
		this.kClusters = clusters;
		this.dimension = dimension;
		this.dataSetSize = dataSetSize;
		this.points = new ArrayList<Point>(dataSetSize);
		this.clustersCentroids = new ArrayList<Point>(kClusters);
		this.clustersMapping = new ArrayList<List<Integer>>(kClusters);
		for (int i = 0; i < kClusters; i++) {
			clustersMapping.add(i, new ArrayList<Integer>());
		}
	}

	// TODO optimize kmeans algorithm for better performance and memory usage
	protected void run() {
		setKRandomPoints();
		// until objective function does not converge... do
		euclideanDistance(points, clustersCentroids);
		// loop while convergence of objective function not satisfied
		// calculate the new centroid points
		//List<Point> newCentroids = new ArrayList<Point>(kClusters);
		int pointsCounter = 0;
		List<Double> sum = new ArrayList<Double>(dimension);
		for(int i = 0; i < dimension; i++){
			sum.add(0d);
		}
		boolean isStopCondition = false;
		while (!isStopCondition) {
			for (int clusterID = 0; clusterID < kClusters; clusterID++) {
				//for (int d = 0; d < dimension; d++) {
					List<Integer> indexList = clustersMapping.get(clusterID);
					for (Integer clusterPoint : indexList) {
						for(int i = 0; i < dimension; i++){
							sum.set(i, sum.get(i) + points.get(clusterPoint).values.get(i));
						}
						//sum += points.get(clusterPoint).values;
						pointsCounter++;
					}
					for(int i = 0; i < dimension; i++){
						sum.set(i, sum.get(i)/pointsCounter);
					}
					Point newCentroid = new Point();
					newCentroid.isVirtual = true;
					newCentroid.values = sum;
					//newCentroid.add(d ,point);
					//newCentroid.add(d, sum / pointsCounter);
					for(int i = 0; i < dimension; i++){
						sum.set(i, 0d);
					}
					pointsCounter = 0;
				//}
				if (!compare(clustersCentroids.get(clusterID), newCentroid)) {
					clustersCentroids.set(clusterID, newCentroid);
					isStopCondition = false;
				} else {
					isStopCondition = true;
				}
				//newCentroid = new ArrayList<Double>(dimension);
			}
			//System.out.println(clustersCentroids.get(0));
			//System.out.println(clustersCentroids.get(1));
			//System.out.println(clustersCentroids.get(2));

			for (int i = 0; i < kClusters; i++) {
				clustersMapping.get(i).clear();
			}

			euclideanDistance(points, clustersCentroids);
		}
	}

	protected boolean compare(Point p1, Point p2) {
		boolean isEqual = true;
		for (int i = 0; i < dimension; i++) {
			isEqual = (p1.values.get(i).doubleValue() == p2.values.get(i).doubleValue()) ? true
					: false;
			if(isEqual == false){
				return isEqual;
			}
		}
		return isEqual;
	}

	/*
	 * Cost function of euclidean distance between points sqrt(sum(xi^2-cj^2))
	 */
	protected void euclideanDistance(List<Point> x, List<Point> c) {
		// TODO calculate euclidean distance between two points in a
		// d-dimensional space
		double sum = 0.0d;
		int indexC = 0;
		double minDistance = 1000d;
		double distance = 0.0d;

		// iterate over all points N
		for (int n = 0; n < x.size(); n++) {
			// iterate over centroids C
			for (int k = 0; k < c.size(); k++) {
				// iterate over dimension
				for (int d = 0; d < dimension; d++) {
					sum += Math.pow(x.get(n).values.get(d)-c.get(k).values.get(d), 2);
							//- Math.pow(c.get(k).values.get(d), 2);
				}
				distance = Math.sqrt(Math.abs(sum));
				/*if (k == 0d) {
					minDistance = distance;
					indexC = k;
				}*/
				if (distance < minDistance) {
					// see after what happens when two points have the same
					// distance
					minDistance = distance;
					indexC = k;
				}
				sum = 0.0d;
			}
			// this is adding duplicate values, check it !!
			clustersMapping.get(indexC).add(n);
			x.get(n).cluster = indexC;
			indexC = 0;
		}
		//this.print();

	}

	/*
	 * TODO Search for methods on picking a good seed
	 */
	protected void setKRandomPoints() {
		Random random = new Random();
		List<Integer> seeds = new ArrayList<Integer>();
		while(seeds.size()<kClusters){
			int ran = random.nextInt(dataSetSize);
			if(!seeds.contains(ran)){
				seeds.add(ran);
			}
		}
		
		for(int k = 0; k < kClusters; k++){
			clustersCentroids.add(k, points.get(seeds.get(k)));
		}
		
	}

	
	protected void readDataSet(File f) {
		/*File f = new File(System.getProperty("user.dir")
				+ "/resources/dataset.txt");*/
		FileReader fr;
		try {
			fr = new FileReader(f);

			BufferedReader br = new BufferedReader(fr);
			String line;
			int skipCount = 0;
			while ((line = br.readLine()) != null) {
				String[] tmps = line.split(",");
				Point point = new Point();
				point.id = Integer.valueOf(tmps[0]);
				//System.out.println(tmps.length);
				for(int i = 1; i < tmps.length; i++){
					//System.out.println(Double.valueOf(tmps[i]));
					point.values.add(i-1,Double.valueOf(tmps[i]));
				}
				point.cluster = -1;
				point.isVirtual = false;
				points.add(skipCount,point);
				skipCount++;
			}

			br.close();
			fr.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	protected void writeResult(File fout) {
		try {
			FileWriter fw = new FileWriter(fout);		
			/*for(int k = 0; k < kClusters; k++){
				List<Integer> indexList = clustersMapping.get(k);
				//Integer[] array = (Integer[]) indexList.toArray();
				//Arrays.sort(array);
				for(int i : indexList){
					fw.write(points.get(i).id+"\t");
				}
				fw.write("\n");
			}*/
			for(int i = 0; i < dataSetSize; i++){
				fw.write(points.get(i).id+"\t"+points.get(i).cluster+"\n");
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		kMeans kmeans = new kMeans(50,100,651);
		File input = new File("D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_CP/kmeansData.csv");
		File output = new File("D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_CP/kmeansResult.txt");
		kmeans.readDataSet(input);
		kmeans.run();
		kmeans.writeResult(output);
	}
}