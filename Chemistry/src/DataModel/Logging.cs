// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


using System;
using Microsoft.Extensions.Logging;

namespace Microsoft.Quantum.Chemistry
{
    public static class Logging
    {
        private static ILoggerFactory loggerFactory = null;

        public static string LogPath
        {
            get => logPath;
            set {
                if (loggerFactory != null)
                {
                    loggerFactory.Dispose();
                    loggerFactory = null;
                }
                else
                {
                    logPath = value;
                }
            }
        }
        private static string logPath;

        public static ILoggerFactory LoggerFactory
        {
            get
            {
                if (loggerFactory == null)
                {
                    loggerFactory = new LoggerFactory()
                        .AddDebug()
                        .AddConsole(true);

                    if (LogPath != null)
                    {
                        System.Console.WriteLine($"Logging to {LogPath}.");
                        loggerFactory
                            .AddFile(LogPath);
                    }
                }

                return loggerFactory;
            }
        }
    }
}
