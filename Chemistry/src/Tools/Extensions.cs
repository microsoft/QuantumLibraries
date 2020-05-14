// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.CommandLine;
using System.CommandLine.Invocation;

namespace Microsoft.Quantum.Chemistry.Tools
{
    public static class Extensions
    {
        public static Command WithHandler(this Command command, ICommandHandler handler)
        {
            command.Handler = handler;
            return command;
        }

        public static Command WithHandler<T1>(this Command command, Action<T1> handler) =>
            command.WithHandler(CommandHandler.Create<T1>(handler));

        public static Command WithHandler<T1, T2>(this Command command, Action<T1, T2> handler) =>
            command.WithHandler(CommandHandler.Create<T1, T2>(handler));
        
        public static Command WithHandler<T1, T2, T3>(this Command command, Action<T1, T2, T3> handler) =>
            command.WithHandler(CommandHandler.Create<T1, T2, T3>(handler));

        public static Command WithHandler<T1, T2, T3, T4>(this Command command, Action<T1, T2, T3, T4> handler) =>
            command.WithHandler(CommandHandler.Create<T1, T2, T3, T4>(handler));

        public static Command WithHandler<T1, T2, T3, T4, T5>(this Command command, Action<T1, T2, T3, T4, T5> handler) =>
            command.WithHandler(CommandHandler.Create<T1, T2, T3, T4, T5>(handler));

        public static Command WithDescription(this Command command, string description)
        {
            command.Description = description;
            return command;
        }

        internal static Broombridge.V0_2.Data WithDefaultMetadata(Broombridge.V0_2.Data data)
        {
            // TODO: make a constant for well-known schema URLs.
            data.Schema = "https://raw.githubusercontent.com/Microsoft/Quantum/master/Chemistry/Schema/broombridge-0.2.schema.json";
            data.Format = Broombridge.V0_2.DefaultFormat;
            data.Generator = new Broombridge.V0_1.Generator
            {
                Source = "qdk-chem",
                Version = typeof(Convert).Assembly.GetName().Version.ToString()
            };
            data.Bibliography = new List<Broombridge.V0_1.BibliographyItem>();
            return data;
        }
    }
}

