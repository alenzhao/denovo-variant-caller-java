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

import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;

/**
 * Abstract Class defining common operations expected of a Caller in this package
 *
 */
public abstract class DenovoCaller {

  /**
   * Execute the DenovoCaller
   *
   * @throws ParseException
   * @throws IOException
   */
  public abstract void execute() throws ParseException, IOException;


  CallHolder parseLine(String line) throws ParseException {
    String[] splitLine = line.split(",");
    if (splitLine.length < 2) {
      throw new ParseException("Could not parse line : " + line, 0);
    }
    String chromosome = splitLine[0];
    Long candidatePosition = Long.valueOf(splitLine[1]);
    CallHolder callHolder = new CallHolder(chromosome, candidatePosition);
    return callHolder;
  }

  /**
   * Write calls to printstream (usually file)
   *
   * @param writer the print stream
   * @param calls A sring describing the calls made, usually csv form
   */
  synchronized void writeCalls(PrintWriter writer, String calls) {
    writer.print(calls);
    writer.flush();
  }

  /**
   * Data class for holding candidate call positions
   */
  class CallHolder {
    final String chromosome;
    final Long position;

    CallHolder(String chromosome, Long position) {
      this.chromosome = chromosome;
      this.position = position;
    }

    @Override
    public String toString() {
      return String.format("<%s,%s>", chromosome, position);
    }
  }
}
