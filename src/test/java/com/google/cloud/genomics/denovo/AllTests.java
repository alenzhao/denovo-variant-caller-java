package com.google.cloud.genomics.denovo;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({NodeTest.class, DenovoBayesNetTest.class, BayesInferTest.class})
public class AllTests {

}
