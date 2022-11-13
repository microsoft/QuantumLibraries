// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using Microsoft.Jupyter.Core;
using Microsoft.Extensions.Logging;
using Microsoft.Quantum.IQSharp.Jupyter;
using System.Linq;
using System;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    public class DisplayableUnitaryOperatorToTextEncoder : IResultEncoder
    {
        private ILogger<DisplayableUnitaryOperatorToTextEncoder> logger;
        public string MimeType => MimeTypes.PlainText;

        public DisplayableUnitaryOperatorToTextEncoder(ILogger<DisplayableUnitaryOperatorToTextEncoder> logger)
        {
            this.logger = logger;
        }

        public EncodedData? Encode(object displayable)
        {
            if (displayable is DisplayableUnitaryOperator op)
            {
                if (op?.Data is null)
                {
                    logger.LogError("Asked to encode a displayable unitary operator, but its data was null. This should not happen.");
                    return null;
                }
                return op.ToString().ToEncodedData();
            }
            else return null;
        }

    }

    public class DisplayableUnitaryOperatorToHtmlEncoder : IResultEncoder
    {
        private readonly ILogger<DisplayableUnitaryOperatorToHtmlEncoder> logger;
        private readonly IConfigurationSource configurationSource;
        public string MimeType => MimeTypes.Html;

        public DisplayableUnitaryOperatorToHtmlEncoder(
            ILogger<DisplayableUnitaryOperatorToHtmlEncoder> logger,
            IConfigurationSource configurationSource
        )
        {
            this.logger = logger;
            this.configurationSource = configurationSource;
        }

        public EncodedData? Encode(object displayable)
        {
            bool IsNearZero(double value) =>
                System.Math.Abs(value) <= configurationSource.TruncationThreshold;

            if (displayable is DisplayableUnitaryOperator op)
            {
                if (op?.Data is null)
                {
                    logger.LogError("Asked to encode a displayable unitary operator, but its data was null. This should not happen.");
                    return null;
                }

                var precision = configurationSource.GetOptionOrDefault(
                    "dump.unitaryPrecision",
                    3
                );

                var rows = new string[op.Data.Dimension];
                var entries = new string[op.Data.Dimension];

                for (var row = 0; row < op.Data.Dimension; ++row)
                {
                    for (var col = 0; col < op.Data.Dimension; ++col)
                    {
                        var format = $"{{0:G{precision}}}";
                        var re = op.Data[row, col, 0];
                        var im = op.Data[row, col, 1];
                        var reFmt = String.Format(format, re);
                        var imFmt = String.Format(format, System.Math.Abs(im)) + "i";
                        if (IsNearZero(re) && IsNearZero(im))
                        {
                            entries[col] = "0";
                        }
                        else if (IsNearZero(im))
                        { 
                            entries[col] = reFmt;
                        }
                        else if (IsNearZero(re))
                        {
                            entries[col] = im < 0.0 ? $"-{imFmt}" : imFmt;
                        }
                        else
                        {
                            entries[col] = $"{reFmt} {(im < 0.0 ? "-" : "+")} {imFmt}";
                        }
                    }

                    rows[row] = String.Join(" & ", entries);
                }

                var outputMatrix = String.Join(" \\\\\n", rows);
                return $@"
                    <table>
                        <tr>
                            <th>Qubit IDs</th>
                            <td>{String.Join(", ", op.Qubits.Select(q => q.Id))}
                        </tr>

                        <tr>
                            <th>Unitary representation</th>
                            <td>$$
                                \left(\begin{{matrix}}
                                    {outputMatrix}
                                \end{{matrix}}\right)
                            $$</td>
                        </tr>
                    </table>
                ".ToEncodedData();
            }
            else return null;
        }

    }

}
