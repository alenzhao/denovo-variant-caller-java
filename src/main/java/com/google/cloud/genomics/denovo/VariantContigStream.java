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
  private long startPosition;
  private long endPosition;
  private List<String> callsetIds;
  private DenovoShared shared;
  
  
  public VariantContigStream(String contig, long startPosition, long endPosition,
      List<String> callsetIds, DenovoShared shared) {
    this.contig = contig;
    this.callsetIds = callsetIds;
    this.startPosition = startPosition;
    this.endPosition = endPosition;
    this.shared = shared;
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
          shared.getDatasetId(),
          null,
          shared.getMaxVariantResults(),
          callsetIds);

    } else if (searchVariantsExecuted.getNextPageToken() != null) {

      requestCount++;
      searchVariantsRequest = DenovoUtil.createSearchVariantsRequest(searchVariantsRequest,
          contig,
          startPosition,
          endPosition,
          shared.getDatasetId(),
          searchVariantsExecuted.getNextPageToken(),
          shared.getMaxVariantResults(),
          callsetIds);
    } else {
      return null;
    }

    if (shared.getDebugLevel() >= 1) {
      System.out.println("Executing Search Variants Request : " + String.valueOf(requestCount));  
    }

    searchVariantsRequestLoaded =
        shared.genomics.variants().search(searchVariantsRequest);
    searchVariantsExecuted = searchVariantsRequestLoaded.execute();
    List<Variant> variants = searchVariantsExecuted.getVariants();
    
    return variants;
  }
}
