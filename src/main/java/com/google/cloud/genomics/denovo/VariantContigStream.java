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

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.Genomics.Variants.Search;
import com.google.api.services.genomics.model.ContigBound;
import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;

import java.io.IOException;
import java.util.List;

/*
 * Creates a Stream of variants for a particular contig
 */
public class VariantContigStream {
  ContigBound contig;
  private int requestCount = 0;
  private SearchVariantsRequest searchVariantsRequest;
  private SearchVariantsResponse searchVariantsExecuted;
  private Search searchVariantsRequestLoaded;
  private String datasetId;
  private Genomics genomics;

  public VariantContigStream(Genomics genomics, ContigBound contig, String datasetId) {
    this.setGenomics(genomics);
    this.contig = contig;
    this.datasetId = datasetId;
    searchVariantsRequest = null;
    searchVariantsRequestLoaded = null;
    searchVariantsExecuted = null;
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
          DenovoUtil.DEFAULT_START_POS,
          contig.getUpperBound(),
          datasetId,
          null);

    } else if (searchVariantsExecuted.getNextPageToken() != null) {

      requestCount++;
      searchVariantsRequest = DenovoUtil.createSearchVariantsRequest(searchVariantsRequest,
          contig,
          DenovoUtil.DEFAULT_START_POS,
          contig.getUpperBound(),
          datasetId,
          searchVariantsExecuted.getNextPageToken());
    } else {
      return null;
    }

    System.out.println("Executing Search Variants Request : " + String.valueOf(requestCount));

    searchVariantsRequestLoaded =
        getGenomics().variants().search(searchVariantsRequest);
    searchVariantsExecuted = searchVariantsRequestLoaded.execute();
    List<Variant> variants = searchVariantsExecuted.getVariants();
    return variants;
  }

  /**
   * @return the genomics
   */
  public Genomics getGenomics() {
    return genomics;
  }

  /**
   * @param genomics the genomics to set
   */
  public void setGenomics(Genomics genomics) {
    this.genomics = genomics;
  }
}
