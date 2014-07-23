package com.google.cloud.genomics.denovo;

import static org.junit.Assert.*;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.api.services.genomics.Genomics;

public class BayesInferTest {

	private static Genomics genomics;
	private static ExperimentRunner expRunner;
	
	@BeforeClass
	public static void setUp() throws Exception {
		String argsString = "stage1 --candidates_file candidate.calls.tmp "+
				"--client_secrets_filename /home/subhodeep/Downloads/client_secrets.json "+
				"--require_all_scopes";
		String[] args = argsString.split(" ");
		
        GenomicsExperiment.cmdLine = new CommandLine();
		GenomicsExperiment.cmdLine.setArgs(args);    	
		genomics = GenomicsExperiment.buildGenomics(GenomicsExperiment.cmdLine).get();

		expRunner = new ExperimentRunner(genomics,GenomicsExperiment.cmdLine);
	}

	@AfterClass
	public static void tearDown() throws Exception {
	}

	@Test 
	public void testGenomicsIsNotNull() {
		assertTrue(genomics != null);
	}

	@Test 
	public void testExpRunnerIsNotNull() {
		assertTrue(expRunner != null);
	}

	
}
