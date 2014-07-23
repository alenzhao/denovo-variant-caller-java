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

import com.google.api.client.auth.oauth2.Credential;
import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp;
import com.google.api.client.extensions.jetty.auth.oauth2.LocalServerReceiver;
import com.google.api.client.googleapis.auth.oauth2.GoogleAuthorizationCodeFlow;
import com.google.api.client.googleapis.auth.oauth2.GoogleClientSecrets;
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport;
import com.google.api.client.http.HttpRequest;
import com.google.api.client.http.HttpRequestInitializer;
import com.google.api.client.http.javanet.NetHttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.jackson2.JacksonFactory;
import com.google.api.client.util.Lists;
import com.google.api.client.util.store.FileDataStoreFactory;
import com.google.api.services.genomics.Genomics;
import com.google.common.base.Optional;

import org.kohsuke.args4j.CmdLineException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.security.GeneralSecurityException;
import java.util.List;


/**
 * Placeholder for running all Genomics Experiments.
 */
public class GenomicsExperiment {

  private static final String APPLICATION_NAME = "Google-GenomicsDeNovoCaller/1.0";
  private static final java.io.File DATA_STORE_DIR =
      new java.io.File(System.getProperty("user.home"), ".store/genomics_denovo_caller");
  private static final String DEVSTORAGE_SCOPE =
      "https://www.googleapis.com/auth/devstorage.read_write";
  private static final String GENOMICS_SCOPE = "https://www.googleapis.com/auth/genomics";
  private static final JsonFactory JSON_FACTORY = JacksonFactory.getDefaultInstance();
  private static final String ROOT_URL = "https://www.googleapis.com/genomics/v1beta";


  private static FileDataStoreFactory dataStoreFactory;
  private static NetHttpTransport httpTransport;
  static CommandLine cmdLine;
  static ExperimentRunner expRunner;

  private static GoogleClientSecrets loadClientSecrets(String clientSecretsFilename) {
    File f = new File(clientSecretsFilename);
    if (f.exists()) {
      try {
        InputStream inputStream = new FileInputStream(new File(clientSecretsFilename));
        return GoogleClientSecrets.load(JSON_FACTORY, new InputStreamReader(inputStream));
      } catch (Exception e) {
        System.err.println("Could not load client_secrets.json");
      }
    } else {
      System.err.println("Client secrets file " + clientSecretsFilename + " does not exist."
          + "  Visit https://developers.google.com/genomics to learn how"
          + " to install a client_secrets.json file.  If you have installed a client_secrets.json"
          + " in a specific location, use --client_secrets_filename <path>/client_secrets.json.");
    }
    return null;
  }

  private static Genomics buildService(final Credential credential) {
    return new Genomics.Builder(httpTransport, JSON_FACTORY, credential)
        .setApplicationName(APPLICATION_NAME)
        .setRootUrl(ROOT_URL)
        .setServicePath("/")
        .setHttpRequestInitializer(new HttpRequestInitializer() {
          @Override
          public void initialize(HttpRequest httpRequest) throws IOException {
            credential.initialize(httpRequest);
            httpRequest.setReadTimeout(60000); // 60 seconds
          }
        })
        .build();
  }

  private static Credential authorize(List<String> scopes) throws Exception {
    GoogleClientSecrets clientSecrets = loadClientSecrets(cmdLine.clientSecretsFilename);
    if (clientSecrets == null) {
      return null;
    }

    GoogleAuthorizationCodeFlow flow = new GoogleAuthorizationCodeFlow.Builder(httpTransport,
        JSON_FACTORY, clientSecrets, scopes).setDataStoreFactory(dataStoreFactory).build();
    return new AuthorizationCodeInstalledApp(flow, new LocalServerReceiver()).authorize(
        "user" + scopes.size());
  }



  public static void main(String[] args) throws IOException {
    System.out.println("-------- Starting Genomics Experiment ---------");

    cmdLine = new CommandLine();

    try {
      // Parse the command line
      cmdLine.setArgs(args);
      Genomics genomics = buildGenomics(cmdLine).get();

      expRunner = new ExperimentRunner(genomics, cmdLine);

      // Entry point for all Experiments
      executeExperiment(cmdLine.stageId);

    } catch (Exception e) {
      cmdLine.printHelp(e.getMessage() + "\n", System.err);
      e.printStackTrace();
      return;
    }
  }

  public static Optional<Genomics> buildGenomics(CommandLine cmdLine) throws CmdLineException,
      GeneralSecurityException, IOException, Exception {

    // Authorization
    List<String> scopes = Lists.newArrayList();
    scopes.add(GENOMICS_SCOPE);
    if (cmdLine.requireAllScopes) {
      scopes.add(DEVSTORAGE_SCOPE);
    }

    httpTransport = GoogleNetHttpTransport.newTrustedTransport();
    dataStoreFactory = new FileDataStoreFactory(DATA_STORE_DIR);
    Credential credential = authorize(scopes);
    if (credential == null) {
      return Optional.absent();
    }

    try {
      credential.refreshToken();
    } catch (NullPointerException e) {
      System.err.append(
          "Couldn't refresh the OAuth token. Are you using a different client secrets file?\n"
          + "If you want to use a different file, first clear your stored credentials: "
          + "http://google-genomics.readthedocs.org/en/latest/api-client-java/resetting_auth.html \n\n");
      return Optional.absent();
    }

    Genomics genomics = buildService(credential);

    return Optional.of(genomics);
  }

  private static void executeExperiment(String stage_id) throws IllegalAccessException,
      IllegalArgumentException, InvocationTargetException, NoSuchMethodException {
    Method[] methods = ExperimentRunner.class.getDeclaredMethods();

    Method methodMatch = null;
    for (Method method : methods) {
      if (method.getName().equals(stage_id)) {
        methodMatch = method;
        break;
      }
    }
    if (methodMatch == null) {
      throw new NoSuchMethodException("No matching method found for Experiment : " + stage_id);
    } else {
      methodMatch.invoke(expRunner);
    }
  }

}
