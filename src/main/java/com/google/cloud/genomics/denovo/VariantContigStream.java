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

import com.google.api.services.genomics.model.SearchVariantsRequest;
import com.google.api.services.genomics.model.SearchVariantsResponse;
import com.google.api.services.genomics.model.Variant;

import java.io.IOException;
import java.math.BigInteger;
import java.util.List;

/**
 * Creates a Stream of variants for a particular contig
 */
public class VariantContigStream {
  private int requestCount = 0;
  private String nextPageToken;
  SearchVariantsRequest request;
  private DenovoShared shared;

  /**
   * @param contig chromosome
   * @param startPosition
   * @param endPosition
   * @param callsetIds list of callset ids which belong to trio
   * @param shared shared project state
   */
  public VariantContigStream(String contig, long startPosition, long endPosition,
      List<String> callsetIds, DenovoShared shared) {
    this.request = new SearchVariantsRequest()
        .setContig(contig)
        .setCallsetIds(callsetIds)
        .setStartPosition(startPosition)
        .setEndPosition(endPosition)
        .setDatasetId(shared.getDatasetId())
        .setMaxResults(BigInteger.valueOf(shared.getMaxVariantResults()));
    this.shared = shared;
  }

  /**
   * @return if Stream has more objects
   */
  public boolean hasMore() {
    return requestCount == 0 || nextPageToken != null;
  }

  /**
   * @return get variants
   * @throws IOException API hangups
   */
  public List<Variant> getVariants() throws IOException {

    requestCount++;
    request.setPageToken(nextPageToken);

    shared.getLogger().finer(
      String.format("Executing Search Variants Request : " + String.valueOf(requestCount)));

    SearchVariantsResponse response = shared.getGenomics().variants().search(request).execute();

    nextPageToken = response.getNextPageToken();
    return response.getVariants();
  }
}