/*
 *Copyright 2014 Google Inc. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */
package com.google.cloud.genomics.denovo;

import static com.google.cloud.genomics.denovo.DenovoUtil.TrioMember.CHILD;

import com.google.api.services.genomics.model.Call;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.Variant;
import com.google.cloud.genomics.denovo.DenovoUtil.Chromosome;
import com.google.cloud.genomics.denovo.VariantsBuffer.PositionCall;
import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.Lists;

import org.javatuples.Pair;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Makes Denovo calls by examining variants at candidate position and checking mendelian 
 * inheritance rules 
 */
public class VariantCaller extends DenovoCaller {

  private final DenovoShared shared;
  private final AtomicInteger variantCounter = new AtomicInteger();
  
  public VariantCaller(DenovoShared shared){
    this.shared = shared;
    
  }
  /* (non-Javadoc)
   * @see com.google.cloud.genomics.denovo.DenovoCaller#execute()
   */
  @Override
  public void execute() throws IOException {
    
    shared.getLogger().info("---- Starting Variant Caller ----");


    // Define Experiment Specific Constant Values
    final File outputFile = DenovoUtil.getNormalizedFile(shared.getOutputFileName());
    shared.getLogger().fine(String.format("Output File : %s", outputFile.getAbsolutePath()));
    
    List<ContigBound> allContigBounds = shared.getGenomics().variants().getSummary()
    .setDatasetId(shared.getDatasetId())
    .setDisableGZipContent(true)
    .execute()
    .getContigBounds();

    
    // Open File Outout handles
    try (PrintWriter callWriter = new PrintWriter(outputFile);) {

      /* Get a list of all the contigs */
      List<ContigBound> contigBounds = FluentIterable
          .from(allContigBounds)
          .filter(new Predicate<ContigBound>() {

            @Override
            public boolean apply(ContigBound cb) {
              return shared.getChromosomes().contains(
                  Chromosome.valueOf(cb.getContig().toUpperCase()));
            }
          }).toList();

      ExecutorService executor = Executors.newFixedThreadPool(shared.getNumThreads());
      /* Iterate through each contig and do variant filtering for each contig */
      for (ContigBound contigBound : contigBounds) {
        Long startContigPos = shared.getStartPosition() == null ? 1L : shared.getStartPosition();
        Long endContigPos = shared.getEndPosition() == null 
            ? contigBound.getUpperBound() : shared.getEndPosition();
        Long strideLength = (endContigPos - startContigPos) / shared.getNumThreads();
        
        shared.getLogger().info("Processing Chromosome : " + contigBound.getContig());
        
        for (int threadIdx = 0 ; threadIdx < shared.getNumThreads() ; threadIdx++) {
          long start = startContigPos + threadIdx * strideLength;
          long end = threadIdx == shared.getNumThreads() - 1 ? endContigPos : start + strideLength - 1;
          Runnable worker = new SimpleDenovoRunnable(callWriter, contigBound.getContig(),
              start, end);
          executor.execute(worker);
        }
      }

      executor.shutdown();
      while (!executor.isTerminated()) {
      }
      shared.getLogger().info("---- Variant caller terminated ----");
    }
  }
  
  /**
   * Run caller through region and record denovo calls if any
   * @param callWriter print stream
   * @param contig chromosome
   * @param startPosition 
   * @param endPosition
   * @throws IOException API hangups
   */
  void callSimpleDenovo(PrintWriter callWriter, String contig, 
      Long startPosition, Long endPosition)
      throws IOException {
    
    // Create new buffer object for storing retreived variants
    VariantsBuffer vbuffer = new VariantsBuffer();
    
    // Create a stream for retreiving variants
    VariantContigStream variantContigStream = new VariantContigStream(contig,
        startPosition, 
        endPosition,
        Lists.newArrayList(shared.getPersonToCallsetIdMap().values()),
        shared);

    // Keep retreiving variants
    while (variantContigStream.hasMore()) {
      StringBuilder builder = new StringBuilder();

      // Get a fresh batch of variants and filter those without calls
      Optional<List<Variant>> variantsFromStream = 
          Optional.fromNullable(variantContigStream.getVariants());
      
      if (!variantsFromStream.isPresent()) {
        return; 
      }
  
      // Filter variants with  missing calls
      Iterable<Variant> variants = FluentIterable
        .from(variantsFromStream.get())
        .filter(new Predicate<Variant>() {
          @Override
          public boolean apply(Variant variant) {
            return variant.getCalls() != null;
          }});

      
      for (Variant variant : variants) {
        
        // Logging related
        synchronized (this) {
          if (variantCounter.get() % DenovoUtil.VARIANT_INFO_STRIDE == 0
              && variantCounter.get() > 0) {
            shared.getLogger().info(
                String.format("%d Variant candidates processed", variantCounter.get()));
          }
          variantCounter.getAndIncrement();
        }
        
        // Push into queue
        for (Call call : variant.getCalls()) {
          vbuffer.checkAndAdd(shared.getCallsetIdToPersonMap().get(call.getCallsetId()), 
              Pair.with(variant, call));
        }
        // Try to process buffer elements eagerly
        while(vbuffer.canProcess()) {
          Optional<PositionCall> nextCall = Optional.fromNullable(vbuffer.retrieveNextCall());
          if (nextCall.isPresent()) {
            if (nextCall.get().isDenovo()) {
              builder.append(String.format("%s,%d,%s%n", contig, nextCall.get().getPosition(),
                  nextCall.get()));

              // Logging
              shared.getLogger().fine(String.format("%s,%d,%s", contig, 
                  nextCall.get().getPosition(), nextCall.get()));
            }
          }
          vbuffer.pop(CHILD);
        }
      }
      writeCalls(callWriter, builder.toString());
    }
    
    // Flush remaining buffer
    StringBuilder builder = new StringBuilder();
    while (!vbuffer.isEmpty(CHILD)) {
      Optional<PositionCall> nextCall = Optional.fromNullable(vbuffer.retrieveNextCall());
      if (nextCall.isPresent()) {
        if (nextCall.get().isDenovo()) {
          builder.append(
              String.format("%s,%d,%s%n", contig, nextCall.get().getPosition(),
                  nextCall.get()));
          
          // Logging
          shared.getLogger().fine(String.format("%s,%d,%s", contig, 
              nextCall.get().getPosition(), nextCall.get()));
        }
      }
      vbuffer.pop(CHILD);
    }
    writeCalls(callWriter, builder.toString());
  }

  /**
   * Wrapper for making multithreaded calls
   */
  private class SimpleDenovoRunnable implements Runnable {

    private final PrintWriter writer;
    private final String contig;
    private final Long startPos;
    private final Long endPos;
    private int numTries = 0;

    public SimpleDenovoRunnable(PrintWriter writer, String contig, 
        Long startPosition, Long endPosition) {
      this.writer = writer;
      this.contig = contig;
      this.startPos = startPosition;
      this.endPos = endPosition;
    }

    @Override
    public void run() {
      try {
        numTries++;
        callSimpleDenovo(writer, contig, startPos, endPos);
      } catch (IOException e) {
        e.printStackTrace();
        if (numTries < shared.getMaxApiRetries()) {
          System.err.printf("Attempt #%d : contig %s%n", numTries + 1, contig);
          run();
        } else {
          System.err.printf("Failed to run contig : %s%n", contig);
        }
      }
    }
  }
}
