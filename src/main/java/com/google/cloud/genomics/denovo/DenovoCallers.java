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

/**
 * Static Factory for Creating DenovoCaller classes 
 */
public class DenovoCallers {

  
  private DenovoCallers() {
    throw new AssertionError("Tried to instantiate non-instantiable class");
  }
  /**
   * create a new VariantCaller which makes denovo calls by looking at the information from variants
   * in VCF files 
   * 
   * @param shared Shared parameters in project
   * @return Caller that uses variants API
   */
  public static DenovoCaller getVariantCaller(DenovoShared shared) {
    // TODO(smoitra): Auto-generated method stub
    return new VariantCaller(shared);
  }

  /**
   * Create a new ReadCaller which makes denovo calls by performing bayesian inference after 
   * collecting all the reads at a particular position
   *  
   * @param shared shared parameters in project
   * @return Caller that uses reads API
   */
  public static DenovoCaller getReadCaller(DenovoShared shared) {
    // TODO(smoitra): Auto-generated method stub
    return new ReadCaller(shared);
  }
}
