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

import static com.google.cloud.genomics.denovo.DenovoUtil.debugLevel;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.Genomics.Variants.Search;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;

import java.io.IOException;
import java.util.List;

/*
 * Creates a Stream of variants for a particular contig
 */
public class VariantContigStream {
  private String contig;
  private int requestCount = 0;
  private SearchVariantsRequest searchVariantsRequest;
  private SearchVariantsResponse searchVariantsExecuted;
  private Search searchVariantsRequestLoaded;
  private String datasetId;
  private Genomics genomics;
  private long maxVariantResults;
  private long startPosition;
  private long endPosition;
  private List<String> callsetIds;
  
  public VariantContigStream(Genomics genomics, String contig, String datasetId, 
      long maxVariantResults, long startPosition, long endPosition, List<String> callsetIds) {
    this.genomics = genomics;
    this.contig = contig;
    this.datasetId = datasetId;
    this.maxVariantResults = maxVariantResults;
    this.startPosition = startPosition;
    this.endPosition = endPosition;
    this.callsetIds = callsetIds;
  }

  public boolean hasMore() {
    if (searchVariantsRequest == null || searchVariantsExecuted.getNextPageToken() != null) {
      return true;
    } else {
      return false;
    }
  }

  public List<Variant> getVariants() throws IOException {

    if (searchVariantsRequest == null) {
      requestCount++;
      searchVariantsRequest = DenovoUtil.createSearchVariantsRequest(null,
          contig,
          startPosition,
          endPosition,
          datasetId,
          null,
          maxVariantResults,
          callsetIds);

    } else if (searchVariantsExecuted.getNextPageToken() != null) {

      requestCount++;
      searchVariantsRequest = DenovoUtil.createSearchVariantsRequest(searchVariantsRequest,
          contig,
          startPosition,
          endPosition,
          datasetId,
          searchVariantsExecuted.getNextPageToken(),
          maxVariantResults,
          callsetIds);
    } else {
      return null;
    }

    if (debugLevel >= 1) {
      System.out.println("Executing Search Variants Request : " + String.valueOf(requestCount));  
    }

    searchVariantsRequestLoaded =
        genomics.variants().search(searchVariantsRequest);
    searchVariantsExecuted = searchVariantsRequestLoaded.execute();
    List<Variant> variants = searchVariantsExecuted.getVariants();
    
    return variants;
  }
}
