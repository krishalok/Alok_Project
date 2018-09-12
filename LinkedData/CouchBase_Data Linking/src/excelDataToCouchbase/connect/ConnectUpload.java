package excelDataToCouchbase.connect;

import java.util.List;

import com.couchbase.client.java.Bucket;
import com.couchbase.client.java.Cluster;
import com.couchbase.client.java.CouchbaseCluster;
import com.couchbase.client.java.document.JsonDocument;
import com.couchbase.client.java.env.CouchbaseEnvironment;
import com.couchbase.client.java.env.DefaultCouchbaseEnvironment;

public class ConnectUpload {
	CouchbaseEnvironment environment = DefaultCouchbaseEnvironment
			.builder()
			.connectTimeout(50000)
			.kvTimeout(25000)
			.queryTimeout(25000)
			.disconnectTimeout(2500000)
			.managementTimeout(25000)
			.socketConnectTimeout(25000)
			.build();
	String ip = "";
	Cluster cluster = CouchbaseCluster.create(environment, ip);
    	Bucket bucket = cluster.openBucket("Bucket_Name","");
    
    	public void connectAndUpload(List<JsonDocument> list) throws Exception {
		list.parallelStream().forEach(doc -> {
			bucket.upsert(doc);
		});
	}
    	public void close(){
    		bucket.close();
    	}
}
